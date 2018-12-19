##### Import the text file #####
beer <- read_delim("Ratebeer.txt",delim = ":", col_names = F)
table(beer$X1)

##### Extract the Features #####
name <- beer[seq(from=1,to=37486579,by=13),]
beerID <- beer[seq(from=2,to=37486579,by=13),]
brewerID <- beer[seq(from=3,to=37486579,by=13),]
ABV <- beer[seq(from=4,to=37486579,by=13),]
style <- beer[seq(from=5,to=37486579,by=13),]
appearance <- beer[seq(from=6,to=37486579,by=13),]
aroma <- beer[seq(from=7,to=37486579,by=13),]
palate <- beer[seq(from=8,to=37486579,by=13),]
taste <- beer[seq(from=9,to=37486579,by=13),]
overall <- beer[seq(from=10,to=37486579,by=13),]
time <- beer[seq(from=11,to=37486579,by=13),]
profileName <- beer[seq(from=12,to=37486579,by=13),]
revText <- beer[seq(from=13,to=37486579,by=13),]

# appearance, aroma, palate, taste, overall: extract numeric

appearance <- as.data.frame(cbind(appearance$X1, str_split_fixed(appearance$X2, pattern = "/", n = 2)))
aroma <- as.data.frame(cbind(aroma$X1, str_split_fixed(aroma$X2, pattern = "/", n = 2)))
palate <- as.data.frame(cbind(palate$X1, str_split_fixed(palate$X2, pattern = "/", n = 2)))
taste <- as.data.frame(cbind(taste$X1, str_split_fixed(taste$X2, pattern = "/", n = 2)))
overall <- as.data.frame(cbind(overall$X1, str_split_fixed(overall$X2, pattern = "/", n = 2)))

table(aroma$V3)
aroma$V3 <- as.numeric(10)

table(appearance$V3)
appearance$V3 <- as.numeric(5)

table(palate$V3)
palate$V3 <- as.numeric(5)

table(taste$V3)
taste$V3 <- as.numeric(10)

table(overall$V3)
overall$V3 <- as.numeric(20)

# Scale the numeric values from 0 to 1
aroma$V4 <- as.numeric(as.character(aroma$V2)) / aroma$V3
appearance$V4 <- as.numeric(as.character(appearance$V2)) / appearance$V3
palate$V4 <- as.numeric(as.character(palate$V2)) / palate$V3
taste$V4 <- as.numeric(as.character(taste$V2)) / taste$V3
overall$V4 <- as.numeric(as.character(overall$V2)) / overall$V3

# Create a final dataframe to manipulate
beers <- data.frame(beerID = beerID$X2, brewerID = brewerID$X2, name = name$X2, ABV = ABV$X2, time = time$X2, appearance = appearance$V4, aroma = aroma$V4, palate = palate$V4, taste = taste$V4, overall = overall$V4, style = style$X2, profileName = profileName$X2, revText = revText$X2)

##### Write the final file to be analyzed #####
write_csv(beers,"beerReviews.csv", col_names = T)
write_csv(x, path, na = "NA", append = FALSE, col_names = !append)
