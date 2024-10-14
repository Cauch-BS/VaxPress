#!/bin/bash

sudo apt-get update && sudo apt-get upgrade -y
sudo apt-get install -y python3-pip python3-venv

# clone vaxpress
git clone https://github.com/Cauch-BS/VaxPress.git
cd VaxPress 
python3 -m venv .venv
source .venv/bin/activate && pip3 install . 
source .venv/bin/activate && pip3 install linearpartition-unofficial linearfold-unofficial