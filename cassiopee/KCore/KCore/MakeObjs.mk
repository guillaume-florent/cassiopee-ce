#!/bin/sh
# ========================================================================
# Enleve les r�pertoires pour ne garder que les noms de fichier
# ========================================================================
#
OBJST=$1

OBJS=`echo ${OBJST} | sed -e "s;[A-Za-z1-90_]*/;;g"`
echo ${OBJS}
