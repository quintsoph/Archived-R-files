#!/usr/bin/env python

import sys

TEDfilename = sys.argv[1]
mutationsfilename = sys.argv[2]
outputfilename = sys.argv[3]

patientlist = []
mutspatlist = []
commonlist = []

with open(TEDfilename, 'r') as TEfile:
  header = TEfile.readline()
  for line in TEfile:
    linesp = line.rstrip('\n').split('\t')
    patient = linesp[0]
    patientlist.append(patient)

with open(mutationsfilename, 'r') as mutfile:
  patient = mutfile.readline().rstrip('\n')
  patient = patient.split('\t')
  mutspatlist.append(patient)

for sublist in mutspatlist:
  for patient in sublist:
    if patient in patientlist:
      commonlist.append(patient)

output = open(outputfilename, 'w+')
for patient in commonlist:
  output.write(patient + '\n')

output.close()
