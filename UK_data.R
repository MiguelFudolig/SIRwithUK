library(tidyverse)
library(XML)
library(splitstackshape)

count.fields("data_2023-Jan-26.csv",sep=",")

uk <- read.csv("data_2023-Jan-26.csv")

uk |> select(variants) |> slice_head(n=1) |> separate_rows()


read.table(text=uk$variants,sep=c(",")) |> read.table(sep=c(":")) 

# uk |>
#   mutate(variants=str_replace_all(variants,"\\*|\\[|\\{|\\]|\\}|\\'","")) |> 
#   mutate(to=strsplit(variants,",|:")) |> 
#   unnest(to) |> 
#   mutate(ind=rep(c(1,2),length.out=n())) |> 
#   group_by(areaType,areaName,areaCode,date,ind) |>
#   mutate(id=row_number()) |> 
#   spread(ind,to) |> 
#   select(-id, -variants) |> 
#   rename("A"="1","B"="2") |> 
#   ungroup() |>
#   mutate(A=str_trim(A)) -> uk_clean
# 
# 
#   
# 
# uk_clean|> 
#   group_by(areaType,areaName,areaCode,date) |> 
#   pivot_wider(names_from=A,values_from = B, values_fn=list) |> 
#   unnest(c(variant,cumWeeklySequenced,newWeeklyPercentage)) |> 
#   mutate(variant=str_trim(variant)) |> 
#   mutate_at(c("cumWeeklySequenced","newWeeklyPercentage"),as.numeric) -> uk_clean

#OR

uk |>
  mutate(variants=str_replace_all(variants,"\\*|\\[|\\{|\\]|\\}|\\'",""))|> 
  separate_longer_delim(cols=variants,delim=",") |> 
  mutate(variants=str_trim(variants)) |> 
  separate(variants,into=c("variables", "values"), sep=":") |> 
  group_by(areaType,areaName,areaCode,date) |> 
  pivot_wider(id_cols= c(areaType,areaName,areaCode,date),names_from="variables",values_from="values", values_fn=list) |> 
  unnest(c(variant,cumWeeklySequenced,newWeeklyPercentage)) |> 
  mutate_at(c("cumWeeklySequenced","newWeeklyPercentage"),as.numeric) |> 
  mutate(variant=str_trim(variant))

uk_clean |> write.csv("UK_cleaned_data.csv")


#reading the csv file
uk_clean <- read.csv("UK_cleaned_data.csv")


#variant surveillance from March 14, 2021 to Dec. 31, 2021
uk_clean |>mutate_at(vars(date),as.Date)|>
  filter(date > "2021-03-14" & date < "2021-12-31") |> 
  filter(!(grepl("22",variant))) |> 
  ggplot(aes(x=date,
             y=newWeeklyPercentage,
             group=variant,
             color=variant,
             fill=variant)) + geom_point(size=2) + geom_line(size=2) + 
  theme_bw() + xlab("Date") + ylab("Percentage of sequenced samples") +
  ggtitle("2021 Sequenced samples") ->plot2021 
  ggsave(plot2021,filename = "2021 Sequenced samples.png", width=12,height=9, units="in")

#variant surveillance from March 14, 2021 to Dec. 31, 2022
uk_clean |>mutate_at(vars(date),as.Date)|>
  filter(date > "2022-01-01" & date < "2022-12-31") |> 
  filter(grepl("22",variant)|grepl("Omicron",variant)) |> 
  ggplot(aes(x=date,
             y=newWeeklyPercentage,
             group=variant,
             color=variant,
             fill=variant)) + geom_point(size=2) + geom_line(size=2) + 
  theme_bw() + xlab("Date") + ylab("Percentage of sequenced samples") +
  ggtitle("2022 Sequenced samples")-> plot2022
  ggsave(plot2022,filename = "2022 Sequenced samples.png", width=12,height=9, units="in")



uk_clean |>mutate_at(vars(date),as.Date) |> 
  filter(grepl("Alpha",variant)|grepl("Delta",variant)) |> 
  filter(date > "2021-03-14" & date < "2021-11-1") |> 
  select(-cumWeeklySequenced) |> 
  group_by(areaType,areaName,areaCode,date) |> 
  pivot_wider(names_from = variant, values_from=newWeeklyPercentage) |> 
  rename("Alpha"="V-20DEC-01 (Alpha)","Delta_B1"="V-21APR-02 (Delta B.1.617.2)") |> 
  mutate(prop = ((Alpha+0.00001)/(Delta_B1+0.00001)), logprop= log(((Alpha+0.00001)/(Delta_B1+0.00001))),totalpropA=Alpha/100,totalpropD=Delta_B1/100) |> 
  ungroup()-> uk_withproportions

#log ratio of the two strains with proportions
uk_withproportions |> ggplot(aes(x=date)) + 
  theme_bw() + 
  geom_point(aes(y=logprop),size=2, color="red") + 
  geom_line(aes(y=logprop),size=2, color="red") +
  geom_point(aes(y=totalpropA*10),size=2) + geom_line(aes(y=totalpropA*10),size=2) + 
  geom_point(aes(y=totalpropD*10),color="blue") + geom_line(aes(y=totalpropD*10),color="blue") + 
  scale_y_continuous(
    name="Log(Proportion)",
    sec.axis = sec_axis(trans=~./10, name="Proportions")
  )


#the proportions at any given time
uk_withproportions |>  ggplot(aes(x=date)) + theme_bw() + geom_point(aes(y=totalpropA),size=2) + geom_line(aes(y=totalpropA),size=2) + 
  geom_point(aes(y=totalpropD),color="blue") + geom_line(aes(y=totalpropD),color="blue")


#with week numbers as x-axis

uk_withproportions |> mutate(week=row_number()-1) |>
  filter(week <=15) |> ggplot(aes(x=week)) + 
  theme_bw() + theme(aspect.ratio=1)+
  geom_point(aes(y=logprop),size=2, color="red") + 
  geom_line(aes(y=logprop),size=2, color="red") +
  geom_point(aes(y=totalpropA*10),size=2) + geom_line(aes(y=totalpropA*10),size=2) + 
  geom_point(aes(y=totalpropD*10),color="blue") + geom_line(aes(y=totalpropD*10),color="blue") + 
  scale_y_continuous(
    name="Log(Proportion)",
    sec.axis = sec_axis(trans=~./10, name="Proportions")
  ) + 
  scale_x_continuous(name="Week Number",
                     breaks = 0:15) + 
  coord_cartesian(xlim=c(0,15),expand=0) ->plotwp

ggsave(plotwp,filename="logproportion and proportion of alpha and delta covid.png",height=6,units="in")


