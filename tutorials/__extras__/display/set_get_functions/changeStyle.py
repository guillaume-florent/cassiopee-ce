#!/usr/bin/python
# coding: utf-8

r"""changeStyle (array)

CPlot.changeStyle: change CPlot display style"""
import Generator as G
import CPlot
import time

a = G.cart((0, 0, 0), (0.1, 0.1, 0.1), (10, 10, 3))
t = [a]
CPlot.display(t, dim=3, mode=0)

for i in range(10):
    CPlot.changeStyle()
    time.sleep(2)

