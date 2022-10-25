#!/bin/bash

HOST_WEB="10.27.54.11"
USER="${1:-`whoami`}"
DESTINATION_DIR="/profils/$USER/models"

# tar --exclude='./.git' --exclude='./share/*' --exclude='*.tgz' --exclude='sync.sh' -zc ./ | ssh "$USER@$HOST_WEB" "rm -rf $DESTINATION_DIR; mkdir -p $DESTINATION_DIR; tar zx -C $DESTINATION_DIR"
tar --exclude='./.git' --exclude='./share/*' --exclude='*.tgz' --exclude='sync.sh' -zc ./ | ssh "$USER@$HOST_WEB" "mkdir -p $DESTINATION_DIR; tar zx -C $DESTINATION_DIR"
