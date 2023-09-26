sudo apt install git -y
sudo install gcc g++ python3 cmake make -y &
git clone https://github.com/jbeder/yaml-cpp.git &
git clone https://github.com/Subhajit-Roy-Partho/oxDNAcpp.git &
git clone https://github.com/Subhajit-Roy-Partho/oxDNAold.git

cd oxDNAold; mkdir build; cd build; cmake ..; make -j$nproc &
cd ~/yaml-cpp; mkdir build; cd build; cmake ..; make -j$nproc; sudo make install;
cd oxDNAold/main/src