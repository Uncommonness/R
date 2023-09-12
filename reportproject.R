rm(list=ls())

library(tseries)
library(TTR)
library(forecast)
library(dplyr)
library(ggplot2)

nep=read.csv("C:/Users/Son/Desktop/TS_팀플/restnep.csv",header = FALSE,sep=",")
nep.ts=ts(data=as.numeric(nep[,2]),frequency = 12,start = c(1999,6))

#1) 시계열 그림
autoplot(nep.ts, main='TS Plot of NEP', theme='bw')
# 평균이 일정하지 않음
# 분산은 일정한 것처럼 보임 (진폭이 일정해보임)
?ts
cor(1:length(nep.ts),log(nep.ts))

# kpss.test(nep.ts) # 귀무가설 기각, 차분 필요
# adf.test(nep.ts) # 귀무가설 기각, 차분 필요 X

autoplot(log(nep.ts), main='TS Plot of log(NEP)', theme='bw') #혼자 실험해봄
# kpss.test(log(nep.ts)) # 귀무가설 기각, 차분 필요
# adf.test(log(nep.ts)) # 귀무가설 기각, 차분 필요 X

## 시간에 따른 분산변화
# Calculate the variance at different time intervals
# variance <- c()
# for (i in 2:length(nep.ts)) {
#   variance <- c(variance, var(nep.ts[1:i]))
# }
# 
# # Create a time vector corresponding to the time points of your data
# time <- 1:length(variance)
# 
# # Plot the change in variance over time
# plot(time, variance, type = "l", xlab = "Time", ylab = "Variance", main = "Change in Variance over Time")
# 


par(mfrow=c(1,2));plot(nep.ts);plot(log(nep.ts)) #분산안정화 (로그변환)
acf(log(nep.ts), plot=T); pacf(log(nep.ts),plot=T)
kpss.test(log(nep.ts))
adf.test(log(nep.ts))

plot(log(nep.ts));plot(log(nep.ts)-decompose(log(nep.ts))$trend) #추세 제거
acf(na.remove(log(nep.ts)-decompose(log(nep.ts))$trend)); pacf(na.remove(log(nep.ts)-decompose(log(nep.ts))$trend))
kpss.test(na.remove(log(nep.ts)-decompose(log(nep.ts))$trend)) # 귀무가설 채택, 차분 필요X
adf.test(na.remove((log(nep.ts)-decompose(log(nep.ts))$trend))) # 귀무가설 기각, 차분 필요X

plot(nep.ts);plot(log(nep.ts)-decompose(log(nep.ts))$trend-decompose(log(nep.ts))$seasonal) #계절 조정
acf(na.remove(log(nep.ts)-decompose(log(nep.ts))$trend-decompose(log(nep.ts))$seasonal), plot=T); pacf(na.remove(log(nep.ts)-decompose(log(nep.ts))$trend-decompose(log(nep.ts))$seasonal),plot=T)
kpss.test(na.remove(log(nep.ts)-decompose(log(nep.ts))$trend-decompose(log(nep.ts))$seasonal)) # 귀무가설 채택, 차분 필요X
adf.test(na.remove(log(nep.ts)-decompose(log(nep.ts))$trend-decompose(log(nep.ts))$seasonal)) # 귀무가설 기각, 차분 필요X

plot(nep.ts);plot(diff(log(nep.ts)-decompose(log(nep.ts))$trend-decompose(log(nep.ts))$seasonal))#차분
acf(na.remove(diff(log(nep.ts)-decompose(log(nep.ts))$trend-decompose(log(nep.ts))$seasonal)), plot=T); pacf(na.remove(diff(log(nep.ts)-decompose(log(nep.ts))$trend-decompose(log(nep.ts))$seasonal)),plot=T)
kpss.test(na.remove(diff(log(nep.ts)-decompose(log(nep.ts))$trend-decompose(log(nep.ts))$seasonal))) #귀무가설 채택, 차분 필요 X
adf.test(na.remove(diff(log(nep.ts)-decompose(log(nep.ts))$trend-decompose(log(nep.ts))$seasonal))) #귀무가설 기각, 차분 필요 X

auto.arima(diff(log(nep.ts)))
auto.arima(nep.ts) %>% autoplot
auto.arima(log(nep.ts)) %>% autoplot
auto.arima(diff(log(nep.ts))) %>% autoplot

checkresiduals(diff(log(nep.ts)))

plot(decompose(log(nep.ts)))
?decompose

#2) 홀트 윈터스 
hw=HoltWinters(log(nep.ts),seasonal = "additive")
plot(hw)
hw2=exp(hw)
hw
write.csv(as.data.frame(tail(cbind(hw$fitted,hw$x),12)),file="hw.csv")

#3) 로그가법모형 - 분해
log(nep.ts) %>% decompose(type="additive") %>%
  autoplot() +
  ggtitle("log decompose with additive model")

dec.gni=decompose(log(nep.ts),type ="additive")
write.csv(as.data.frame(tail(cbind(dec.gni$trend,dec.gni$seasonal,dec.gni$random,dec.gni$x),18)),file="dec.csv")

#4) 회귀분석
tt=1:length(nep.ts)
reg=lm(log(nep.ts)~tt+I(tt^2)+I(tt^3))
autoplot(nep.ts,main="reg with trend")
par(new=T)
plot(tt,exp(reg$fitted.values),axes=F,ann=F,col="red",type="l")
summary(reg)
write.csv(as.data.frame(tail(cbind(exp(reg$fitted.values),exp(reg$model[,1])),12)),file="reg.csv")

#4-2) 계절성 추가
Q1=(tt%%12)==1
Q2=(tt%%12)==2
Q3=(tt%%12)==3
Q4=(tt%%12)==4
Q5=(tt%%12)==5
Q6=(tt%%12)==6
Q7=(tt%%12)==7
Q8=(tt%%12)==8
Q9=(tt%%12)==9
Q10=(tt%%12)==10
Q11=(tt%%12)==11
regs=lm(log(nep.ts)~tt+I(tt^2)+I(tt^3)+Q1+Q2+Q3+Q4+Q5+Q6+Q7+Q8+Q9+Q10+Q11)
summary(regs)
autoplot(nep.ts,main="reg with smoothing and seasonal")
par(new=T)
plot(tt,exp(regs$fitted.values),axes=F,ann=F,col="red",type="l")
write.csv(as.data.frame(tail(cbind(exp(regs$fitted.values),exp(regs$model[,1])),12)),file="regs.csv")

#5 모형 적합 

## ppt순서
# 분산안정화 - 추세 조정 - 계절조정 (계절아리마로 대체) : 
# 비정상인지 정상인지 판단 -> 차분필요하면 차분 // = acf / pacf 그림 파악

acf(log(nep.ts))
pacf(log(nep.ts))
# => 비정상 시계열 데이터로 보임

par(mfrow=c(1,2))
acf(diff(log(nep.ts)))
pacf(diff(log(nep.ts)))
# => 빠르게 0으로 수렴하지 않는다
# => 주기마다 큰 값을 갖는 것으로 보아 계절arima로 보임

# 지금까지 모형의 식별
# 계절 차분과 차분을 통해 비계절성과 계절성이 있는 arima모형으로 보인다. => 단위근검정필요.
# 계절 차분=>
diff(diff(log(nep.ts)),lag=12)
# 단위근 + kpss, 시계열그림 = 잠정적 모형 선택 /
acf(diff(log(nep.ts),lag=12))
pacf(diff(log(nep.ts),lag=12))  
adf.test(diff(diff(log(nep.ts)),lag=12))
kpss.test(diff(diff(log(nep.ts)),lag=12))
# 모수 추정 방법 선택 : css-ml
X=diff(diff(log(nep.ts)),lag=12)
arima(X,order=c(0,0,0),seasonal = c(2,0,2))
# 모형진단 : 잔차가 WN인가// 과대적합

# 예측 : forecast

# 결국 우리가 사용해야할 모형은 non-seasonal X seasonal 승법모형으로 보임임

autoplot(log(nep.ts), main='TS Plot of NEP', theme='bw')
diff(log(nep.ts),lag=12)
autoplot(diff(log(nep.ts),lag=12))
kpss.test(diff(log(nep.ts),lag=12),null = "Level",lshort = T)
auto.arima(log(nep.ts))
acf(diff(diff(log(nep.ts)),lag=12))
pacf(diff(diff(log(nep.ts)),lag=12))


autoplot(diff(diff(log(nep.ts)),lag=12))
