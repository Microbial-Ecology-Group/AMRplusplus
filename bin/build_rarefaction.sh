git clone https://github.com/cdeanj/rarefactionanalyzer.git
cd rarefactionanalyzer
make
chmod 777 rarefaction
mv rarefaction {$baseDir}/bin/rarefaction
cd ../
rm -rf rarefactionanalyzer
