library(tidyverse)
library(readxl)
library(glue)

# 상권영역
# -----------------------------------------------------------------------------------------------------------
# ID    7자리
# 기준  201810

AREA <- read_csv('상권영역.csv')
colnames(AREA) <- str_remove_all(colnames(AREA), '_')
AREA <- AREA %>% filter(상권구분코드 == 'A')
AREA <- AREA %>% select(-상권구분코드,-상권구분코드명,-형태정보)
AREA <- AREA %>% arrange(시군구코드,행정동코드,상권코드명)
colnames(AREA) <- c('CODE','NAME','LOC_X','LOC_Y','GU','DONG','YEARMONTH')

# 상권코드 비교
code0 <- read_delim('TBSM_TRDAR_SELNG.txt', delim = '|')
code0 <- code0 %>% select(2, 3) %>% distinct()
colnames(code0) <- c('CODE0','NAME')
code0 <- code0 %>% arrange(NAME)
View(code0)

code1 <- AREA  %>% select(CODE, NAME) %>% distinct()
code1 <- code1 %>% arrange(NAME)
colnames(code1) <- c('CODE1','NAME')

CODEBOOK <- code1 %>% left_join(code0)


# 추정매출 
# -----------------------------------------------------------------------------------------------------------
# ID    7자리
# 기준  201804(~201401)

SELNG <- read_csv('상권-추정매출.csv')
colnames(SELNG) <- str_remove_all(colnames(SELNG), '_')
SELNG <- SELNG %>% filter(상권구분코드 == 'A')
SELNG <- SELNG %>% select(-상권구분코드, -상권구분코드명)

SELNG <- SELNG %>% select(1:8, 10, 20:23, 25, 26:28, 78)
# 기준년코드,기준분기코드,상권코드,상권코드명,서비스업종코드,서비스업종코드명,
# 당월매출금액,당월매출건수,주말매출비율,여성매출비율,연령대20매출비율,연령대30매출비율
# 시간대11~14매출비율,시간대14~17매출비율,시간대17~21매출비율,점포수

serviceBook <- SELNG %>% 
  select(서비스업종코드, 서비스업종코드명) %>% distinct() %>% 
  arrange(서비스업종코드명)
serviceBook <- serviceBook %>% 
  filter(서비스업종코드명 %in% 
                   c('양식음식점','일식음식점','중식음식점','한식음식점',
                     '분식전문점','치킨전문점','패스트푸드점',
                     '제과점','커피·음료','호프·간이주점'))

SELNG <- SELNG %>% 
  filter(서비스업종코드명 %in% serviceBook$서비스업종코드명) %>%
  filter(점포수 > 0) 

SELNG <- SELNG %>% filter(기준년코드 >= 2018)
colnames(SELNG) <- c('YEAR','QUAR','CODE','NAME','SVCODE','SVNAME','SALES_MONTH_AMT','SALES_MONTH_CNT',
                     'SALES_PER_WKND', 'SALES_PER_1114', 'SALES_PER_1417', 'SALES_PER_1721','SALES_PER_2124',
                     'SALES_PER_WOMAN', 'SALES_PER_AGE10','SALES_PER_AGE20','SALES_PER_AGE30','N_STORE')


# 소득소비 
# -----------------------------------------------------------------------------------------------------------
# ID    6자리
# 기준  201810(~201401)  

INCN <- read_csv('상권-소득소비-S.csv')
colnames(INCN) <- colnames(INCN) %>% str_remove_all('_') 

INCN <- INCN %>% 
  select(기준년월코드,상권코드,월평균소득금액,지출총금액)

INCN <- INCN %>% filter(기준년월코드 >= 201801)
INCN <- INCN %>% separate(기준년월코드, into = c('YEAR','MONTH'), 4)
colnames(INCN)[3:5] <- c('CODE0','INCM_MONTH_AVG','CNSM_TOT')

# 분기별 요약
INCN_QT <- INCN %>% mutate(MONTH = as.integer(MONTH)) %>% 
  mutate(QUARTER = ifelse(MONTH <= 3, 1, ifelse(MONTH <= 6, 2, ifelse(MONTH <= 9, 3, 4)))) %>%
  group_by(YEAR, QUARTER, CODE0) %>% 
  summarise(IN_MONTH_AVG = mean(INCM_MONTH_AVG), CN_TOT = mean(CNSM_TOT))

INCN_QT <- INCN_QT %>% filter(YEAR == 2018)
INCN_QT$CODE0 <- str_pad(INCN_QT$CODE0, 6, 'left', 0)


# 유동인구
# -----------------------------------------------------------------------------------------------------------
# ID    6자리
# 기준  201810(~201809)   

FLPOP <- read_csv('상권-추정유동인구-S.csv')
colnames(FLPOP) <- str_remove_all(colnames(FLPOP), '_')


# 프로파일링 데이터 엮기(~201801)
load.FLPOP <- function(YY, MM) {
  MM <- str_pad(MM, 2, 'left', 0)
  fname <- glue('원데이터/{YY}년{MM}월/TBSM_TRDAR_FLPOP.txt')
  rawdat <- read_delim(fname, delim = '|') %>% select(1:4,6)
  return(rawdat)
}

rawFLPOP <- NULL
for (m in 1:8) {
  tmp <- load.FLPOP(18, m)
  rawFLPOP <- rbind(rawFLPOP, tmp)
}

colnames(FLPOP) <- c('YM','CODE0','NAME','FLPOP_TOT','FLPOP_WOMAN')
colnames(rawFLPOP) <- c('YM','CODE0','NAME','FLPOP_TOT','FLPOP_WOMAN')
FLPOP <- FLPOP %>% rbind(rawFLPOP)

# 분기별 요약
FLPOP_QT <- FLPOP %>% separate(YM, into = c('YEAR','MONTH'), 4) %>%
  mutate(MONTH = as.integer(MONTH)) %>%
  mutate(QUARTER = ifelse(MONTH <= 3, 1, ifelse(MONTH <= 6, 2, ifelse(MONTH <= 9, 3, 4)))) %>%
  group_by(YEAR, QUARTER, CODE0) %>% 
  summarise(FLPOP_TOT_AVG = mean(FLPOP_TOT), FLPOP_WOMAN_AVG = mean(FLPOP_WOMAN))

FLPOP_QT <- select(FLPOP_QT, -FLPOP_WOMAN_AVG)


# 상주인구 
# ----------------------------------------------------------------------------------------------------------
# ID    6자리
# 기준  201810(~201401) 

REPOP <- read_csv('상권-상주인구-S.csv') 
colnames(REPOP) <- str_remove_all(colnames(REPOP), ' ')
REPOP <- REPOP %>% select(기준년월코드, 상권코드, 상권코드명, 총상주인구수)
colnames(REPOP) <- c('YM','CODE0','NAME','REPOP_TOT')
REPOP <- REPOP %>% filter(YM >= 201801)

# 분기별 요약
REPOP_QT <- REPOP %>% separate(YM, into = c('YEAR','MONTH'), 4) %>%
  mutate(MONTH = as.integer(MONTH)) %>%
  mutate(QUARTER = ifelse(MONTH <= 3, 1, ifelse(MONTH <= 6, 2, ifelse(MONTH <= 9, 3, 4)))) %>%
  group_by(YEAR, QUARTER, CODE0, NAME) %>% 
  summarise(REPOP_TOT_AVG = mean(REPOP_TOT))
REPOP_QT$CODE0 <- REPOP_QT$CODE0 %>% str_pad(6, 'left', 0)

# 인구/소득소비정보 MERGE
POPN <- REPOP_QT %>% left_join(FLPOP_QT) 
POPN <- POPN %>% left_join(INCN_QT)

SELNG <- SELNG %>% rename(QUARTER = QUAR)
POPN$YEAR <- as.integer(POPN$YEAR)


# 최종 데이터셋 *** 
MYDAT <- SELNG %>% left_join(POPN)
MYDAT <- MYDAT %>% na.omit() 

MYDAT <- MYDAT %>% select(YEAR, QUARTER, CODE, CODE0, NAME, SVCODE, SVNAME, starts_with('SALES'), N_STORE,
                          REPOP_TOT_AVG, FLPOP_TOT_AVG, IN_MONTH_AVG, CN_TOT_AVG = CN_TOT)
write.csv(MYDAT, 'DATA2018-FINAL-FOOD.csv', row.names = F, fileEncoding = 'EUC-KR')


