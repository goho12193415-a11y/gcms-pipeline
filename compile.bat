copy .\App2.java .\src\main\java\ru\ac\phyche\gcmsburyak\retentionprediction4\App2.java
copy .\TrainPolar.java .\src\main\java\ru\ac\phyche\gcmsburyak\retentionprediction4\TrainPolar.java
copy .\SecondLevelModelForPolar.java .\src\main\java\ru\ac\phyche\gcmsburyak\retentionprediction4\SecondLevelModelForPolar.java
mkdir models_polar
copy .\db624.svr .\models_polar\db624.svr
copy .\db17.svr .\models_polar\db17.svr
copy .\models\descriptors_info.txt .\models_polar\descriptors_info.txt
call mvn clean compile test
call mvn package
copy .\target\retentionprediction4-0.0.6-jar-with-dependencies.jar .\retentionprediction4-0.0.6-jar-with-dependencies.jar
call java -Xmx1500M -cp retentionprediction4-0.0.6-jar-with-dependencies.jar ru.ac.phyche.gcmsburyak.retentionprediction4.TrainPolar txt2nn .\mlp.txt .\cnn.txt .\models_polar\mlp.nn .\models_polar\cnn.nn
call java -Xmx1500M -cp retentionprediction4-0.0.6-jar-with-dependencies.jar ru.ac.phyche.gcmsburyak.retentionprediction4.TrainPolar txt2nn .\mlpPolar.txt .\cnnPolar.txt .\models_polar\mlpPolar.nn .\models_polar\cnnPolar.nn
call lazbuild retentionindexprediction.lpi


