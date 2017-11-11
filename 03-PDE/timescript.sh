#!/bin/bash

HOSTNAME=$(hostname -s)
TIMESTAMP=$(date +"%Y-%m-%d %H-%M-%S.%N")

echo "$HOSTNAME: $TIMESTAMP"
