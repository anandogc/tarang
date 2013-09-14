#!/bin/bash

find ../../.. -type f  \( -iname "*.h" -or -iname "*.cc" -or -iname "*.txt" \) > /tmp/tarang_source_list
awk -f lc.awk $(cat /tmp/tarang_source_list)
