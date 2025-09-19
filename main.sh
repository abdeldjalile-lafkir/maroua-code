#!/bin/bash
# Send the last commit to the VM after each commit
VM_USER=lafkir_abdeldjalile
VM_HOST=your.vm.ip
VM_PATH=/home/lafkir_abdeldjalile/code

git push origin main  # push commit to remote repo
ssh $VM_USER@$VM_HOST "cd $VM_PATH && git pull origin main"
