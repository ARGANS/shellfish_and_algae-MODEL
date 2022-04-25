#!/bin/bash

HOST_WEB="10.27.54.11"
USER="${1:-`whoami`}"
DESTINATION="/profils/$USER/prj1"

tar --exclude='./.git' --exclude='./share/*' --exclude='*.tgz' --exclude='sync.sh' -zc ./ | ssh "$USER@$HOST_WEB" "rm -rf $DESTINATION; mkdir -p $DESTINATION; tar zx -C $DESTINATION"
