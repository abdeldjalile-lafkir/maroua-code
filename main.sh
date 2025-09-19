#!/bin/bash
# Send the last commit to the VM after each commit
VM_USER=lafkir_abdeldjalile
VM_HOST=34.16.88.151
VM_PATH=/home/lafkir_abdeldjalile/code

git add .
git commit -m "Auto commit from local machine"
git push origin master




ssh $VM_USER@$VM_HOST "cd $VM_PATH && git pull origin master
