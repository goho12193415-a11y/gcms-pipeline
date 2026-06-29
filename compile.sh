cp ./App2.java ./src/main/java/ru/ac/phyche/gcmsburyak/retentionprediction4/App2.java
cp ./TrainPolar.java ./src/main/java/ru/ac/phyche/gcmsburyak/retentionprediction4/TrainPolar.java
cp ./SecondLevelModelForPolar.java ./src/main/java/ru/ac/phyche/gcmsburyak/retentionprediction4/SecondLevelModelForPolar.java
mkdir models_polar
cp ./db624.svr ./models_polar/db624.svr
cp ./db17.svr ./models_polar/db17.svr
cp ./models/descriptors_info.txt ./models_polar/descriptors_info.txt
mvn clean compile test
mvn package
cp ./target/retentionprediction4-0.0.6-jar-with-dependencies.jar ./retentionprediction4-0.0.6-jar-with-dependencies.jar
java -Xmx1500M -cp retentionprediction4-0.0.6-jar-with-dependencies.jar ru.ac.phyche.gcmsburyak.retentionprediction4.TrainPolar txt2nn ./mlp.txt ./cnn.txt ./models_polar/mlp.nn ./models_polar/cnn.nn
java -Xmx1500M -cp retentionprediction4-0.0.6-jar-with-dependencies.jar ru.ac.phyche.gcmsburyak.retentionprediction4.TrainPolar txt2nn ./mlpPolar.txt ./cnnPolar.txt ./models_polar/mlpPolar.nn ./models_polar/cnnPolar.nn
lazbuild retentionindexprediction.lpi


