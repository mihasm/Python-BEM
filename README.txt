README

How to create executable, you may ask? Without further ado, copy and paste this line into your command console:

pyinstaller .\main.py --onefile --paths C:\\Users\\miha\\AppData\\Local\\Programs\\Python\\Python36\\lib\\site-packages\\scipy\\extra-dll --hidden-import='scipy._lib.messagestream' --windowed