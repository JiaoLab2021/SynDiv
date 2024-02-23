#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import gzip

# opener
class MyOpener():
    def __init__(self, filePath):
        # genome
        self.filePath = filePath

        # open file
        self.fileIO = self.open_file()

    # build index
    def read_line(self, line):
        # remove
        line[0] = None
        # read
        lineValue = self.fileIO.readline().strip()
        if lineValue:
            line[0] = lineValue
            return True
        else:
            self.close_file()
            return None

    # open file
    def open_file(self):
        if self.filePath.endswith('.gz'):
            return gzip.open(self.filePath, 'rt', encoding='UTF-8')
        else:
            return open(self.filePath, 'r', encoding='UTF-8')

    # close file
    def close_file(self):
        if not self.fileIO.closed:
            self.fileIO.close()