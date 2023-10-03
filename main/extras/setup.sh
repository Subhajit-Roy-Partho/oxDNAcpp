sudo apt install git tmux -y
sudo apt install gcc g++ python3 cmake make libeigen3-dev htop -y
git clone https://github.com/jbeder/yaml-cpp.git &
git clone https://github.com/Subhajit-Roy-Partho/oxDNAcpp.git &
git clone https://github.com/Subhajit-Roy-Partho/oxDNAold.git &

cd oxDNAold; mkdir build; cd build; cmake ..; make -j$nproc; make romano -j$nproc;cd bin;echo "export PATH=$PATH:$(readlink -f .)" >> ~/.bashrc; echo "shopt -s autocd" >> ~/.bashrc;source ~/.bashrc;
cd ~/yaml-cpp; mkdir build; cd build; cmake ..; make -j$nproc; sudo make install;
cd oxDNAcpp/main/src;mkdir build; cd build; cmake ..; make -j$nproc
