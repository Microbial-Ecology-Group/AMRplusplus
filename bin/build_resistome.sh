git clone https://github.com/cdeanj/resistomeanalyzer.git
cd resistomeanalyzer
make
chmod 777 resistome
mv resistome {$baseDir}/bin/resistome
cd ../
rm -rf resistomeanalyzer
