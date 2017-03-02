#python with functions for minipacket data extraction
import matplotlib.pyplot as plt
import os.path
from numpy import array, zeros, where, arange, full, histogram, mean
import pyfits
from datetime import datetime


files = ['fe55_pwm_45_sdd0_0p08V_sdd1_1p0V.bin']
fits_name = 'test.fits'

#read in data packets from file
def read_into_buffer(filename):
  data = bytearray()
  for i in range(len(filename)):
    buffer = bytearray(os.path.getsize(filename[i]))
    with open(filename[i], 'rb') as f:
      f.readinto(buffer)
      data += bytearray(buffer)
  return data

def init_crcccitt_tab(): #function from lib_crc.c to generate CRC validation table
  P_CCITT = 0x1021 #CRC polynomial
  crc_tabccitt = zeros(256, dtype=int)
  for i in range(256):
    crc = 0x0000
    c = ((i & 0xFFFF) << 8) & 0xFFFF
    for j in range(8):
      if (((crc ^ c) & 0x8000)):
        crc = ((crc << 1) ^ P_CCITT) & 0xFFFF
      else: crc = (crc << 1) & 0xFFFF
      c = (c << 1) & 0xFFFF
      crc_tabccitt[i] = (crc) & 0xFFFF
  return crc_tabccitt

def update_crc_ccitt(table, crc, c): #function from lib_crc.c to update CRC
  short_c = (c & 0xFF)
  #crc = (crc & 0xFFFF)
  tmp = ((crc >> 8) ^ short_c)
  crc = ((crc << 8) ^ table[tmp])
  #print format (c, '02X'), format(short_c, '06X'), format(crc, '06X'), format(tmp, '02X')
  return crc

def crc_check(table, packet, good_crc):
  test_crc = 0x0000
  for i in range(len(packet)):
    #packet is transport frame starting from sequence monotonic to end of packet
    test_crc = (update_crc_ccitt(table, test_crc, packet[i])) & 0xFFFF
    #print format(test_crc, '04X'), format(good_crc, '04X')
  if (test_crc == good_crc):
    return True
  else:
    return False

def UTIL_event_time(mp_time_in, wave_time_in, epoch): 
#calculate time minipacket was created from truncated time
  MAX_DELTA = int(64)
  pkt_time = wave_time_in & 0xFFFFFFFE #strip time quality bit
  mp_time = int(mp_time_in)
  mp_time = int(mp_time/40) #10 bit second in minipacket
  delta = mp_time - (pkt_time & 0x000003FF) #32 bit SCLK
  if (delta > MAX_DELTA):
    delta -= 0x0400
  return pkt_time + delta + epoch #return minipacket time

def telemetry(f):
  l = len(f) #number of bytes in file
  science = bytearray() #create empty byte array to hold LRS minipackets
  hsk = bytearray() #create empty byte array to hold HSK minipackets
  i = 0 #pointer index
  n = 0 #keep track of number of teletry packets
  s = 0 #keep track of number of science packets
  m = 0 #keep track of number of LRS minipackets
  h = 0 #keep track of number of housekeeping packets
  c = 0 #keep track of bad CRC's
  x = 0 #keep track of number of X-ray events
  crc_tabccitt = init_crcccitt_tab() #create CRC table
  while (i < l): #go through file and pull out transport packets looking for sync pattern
    if ((f[i] == 0xFA) and (f[i+1] == 0xF3) and (f[i+2] == 0x34) and (f[i+3] == 0x03)):
      #good sync pattern do stuff
      #crc is bytes 4, 5 indexing from 0
      crc = ((f[i+4]&0xFF)<<8) | ((f[i+5]&0xFF)<<0)
      #CRC validation, if bad move 1 ahead
      #record length, bytes 14, 15 length from first minipacket byte to end of fill data
      length = ((f[i+14]&0xFF)<<8) | ((f[i+15]&0xFF) <<0)
      if (crc_check(crc_tabccitt, f[i+6:i+16+length], crc) == True): #do a bunch of other stuff
        #sequence monotonic is bytes 6, 7
        #first two bits of byte 6 are the source, last 14 are frame count
        source = ((f[i+6]&0xC0)>>6)&0x03
        #frame = ((f[i+6]&0x3F)<<8) | ((f[i+7]&0xFF)<<0)
        #system time is bytes 8, 9, 10, 11 LSB of 11 is quality flag
        sys_time = ((f[i+8]&0xFF)<<24) | ((f[i+9]&0xFF)<<16) | ((f[i+10]&0xFF)<<8) | ((f[i+11]&0xFE))
        #print sys_time
        #check the data quality bit
        if ((f[i+11]&0x01) == 1): #got a bad time flag, need to do stuff in minipacket
          bad_time = True
          print 'Bad Time!'
          #do stuff if bad time flag
        #epoch is bytes 12, 13
        #epoch = ((f[i+12]&0xFF)<<8) | ((f[i+13]&0xFF)<<0)
        if (source == 0x02): #packet is science telemetry
          """ Start Mini Packet Extraction Here"""
          #start mini packet work, create index j which is at first byte of first minipacket
          j = i+16 #minipacket bytes are now relative to j
          while (j < (i + 16 + length)): #only look for minipackets in this frame
            #byte 0 4 MSB are the data type
            mini_type = ((f[j]&0xF0)>>4) #&0x0F
            #byte 0 bottom 4 bits and byte 1 are length
            mini_length = (((f[j]&0x0F)<<8) | ((f[j+1])&0xFF))
            #print mini_type
            if (mini_type == 0x00): #0 means fill data, move to next frame
              break #kill while loop
            elif (mini_type == 2): # HaloSat X-ray packet, continue with packet extraction
              #UPDATE ME
              #bytes 2, 3 are truncated time
              trunc_time = ((f[j+2]&0xFF)<<8) | (f[j+3]&0xFF)
              #byte 4, 5 are segmentation and MSF
              MSF = (f[j+4] >> 5) & 0x03
              if (MSF > 0): print 'MSF = ', MSF
              #if (MSF == 0x00):
              #mini_length += 9
              #elif (MSF == 0x01):
              #  mini_length += 11
              #elif (MSF == 0x02):
              #  mini_length += 13
              #elif (MSF == 0x03):
              #  mini_length += 17
              #print (mini_length - 6)/4
              #print mini_length
              #calculate system time minipacket was created
              mini_time = UTIL_event_time(trunc_time, sys_time, 0)
              #print 'hello(15), sys(32), mini(32)', trunc_time/40, sys_time, mini_time
              #build up new header for mini packet
              #add recalculated system time to minipacket header
              #print ((mini_time >> 24)&0xFF) | ((mini_time >> 16)&0xFF) | ((mini_time >>8)&0xFF) | ((mini_time >>0)&0xFF)
              for k in range(j+6, j+mini_length, 4) :
                #trunc_event_time = ((f[k]&0xFF)<<8) | ((f[k+1]&0xFE)) #15 bit truncated event time
                trunc_event_time = ((f[k]&0x7F)<<8) | ((f[k+1]&0xFE)) #15 bit truncated event time
                #event_time = trunc_event_time | sys_time&0xFFF8000  #UTIL_event_time(trunc_event_time, mini_time, 0)
                event_time = trunc_event_time # | ((sys_time&0xFFFFC00)*40)  #UTIL_event_time(trunc_event_time, mini_time, 0)
                etime.append(event_time)
                event_status = ((f[k+1]&0x01)<<2)+((f[k+2]&0xC0)>>6)
                #print event_status
                stat.append(event_status)
                #stat.append(MSF)
                # 1=X-ray, 3=reset, 4=test
                pulse_height = ((f[k+2]&0x3F)<<8) | ((f[k+3]&0xFF))
                ph.append(pulse_height)
                # print trunc_event_time, pulse_height
            #   for z in range(1, 5):
            #     time_byte = (mini_time >> (32-(8*z))) & 0xFF
            #     header.append(time_byte)
            #     science += header + f[j:j+mini_length]
                #move pointer ahead by length of mini packet
              j += mini_length
              x += (mini_length - 6)/4
              m += 1
            else :
              j += mini_length
          s += 1
          i += 16 + length
          n += 1
        elif (source == 0x03): #packet is Housekeeping telemetry
          hsk += f[i: i+16+length] #store housekeeping packet
          n += 1
          h += 1
          i += 16 + length
      else: #move pointer ahead by 1 if crc is not valid
        i += 1
        c += 1
        #print 'Bad CRC!'
    else: #move on to next byte, sync pattern not found
      #print 'No sync pattern detected! check code/packets!'
      i+= 1

  if (c > 0):
      print 'Bad CRC!'
  print 'Good Transport Frames: ', n
  print 'Bad CRCs: ', c
  print 'LRS packets: ', s
  print 'Housekeeping Packets: ', h
  print 'Science Minipackets', m
  print 'Number of events: ', x
  print 'Total Number of Packets (sanity check): ', s+h
  return science, hsk, m, x #function to search through binary file for telemetry

def mini_extract(f, num_mini, num_x): #function to read in events from minipackets
  l = len(f)
  #print 'Total science minipacket length: ', l
  #print 'Number of bytes left after subtraction: ', l - (6*num_mini) - (4*num_x) - (4*num_mini)
  i = 0 #pointer index
  p = 0 #array storage index
  num_events = num_x
  #num_events = (l - (6 + 4)*num_mini)/4 #total number of X-ray events in minipackets
  #print 'Number of X-ray events: ', num_events
  #create arrays to store event data
  event_time, test, reset = zeros(num_events), zeros(num_events), zeros(num_events)
  X_ray, pulse_height =  zeros(num_events), zeros(num_events)
  while (i < l):
    #first four bytes of a given minipacket are the reconstructed time
    mini_time = ((f[i]&0xFF)<<24) | ((f[i+1]&0xFF)<<16) | ((f[i+2]&0xFF)<<8) | ((f[i+3]&0xFF))
    mini_type = ((f[i+4]&0xF0)>>4)&0x0F
    #print 'LRS', mini_time
    mini_length = ((f[i+4]&0x0F)<<8) | ((f[i+5])&0xFF)
    MSF = (f[i+4+4] >> 5) & 0x03
    if (MSF == 0x00):
      mini_length += 9
    elif (MSF == 0x01):
      mini_length += 11
    elif (MSF == 0x02):
      mini_length += 13
    elif (MSF == 0x03):
      mini_length += 17
    #print mini_length
    if (mini_type == 0x02):
      j = i + 4 + (8 - 2) #create pointer at first X-ray event
      while (j < i + 4 + mini_length):
        trunc_event_time = ((f[j]&0xFF)<<8) | ((f[j+1]&0xFE)) #15 bit truncated event time
        event_time[p] = UTIL_event_time(trunc_event_time, mini_time, 0)
        #print 'tunc_event(15), mini time(32), event(32), type', trunc_event_time/40, mini_time, event_time[p], mini_type
        #check status bits for event
        if ((f[j+1]&0x01) == 0x01): #event is from a test input
          test[p] = 1
        elif ((((f[j+2]&0xC0)>>6) & 0x03) == 0x03): #event is reset pulse
          reset[p] = 1
        elif ((((f[j+2]&0xC0)>>6) & 0x01) == 0x01): #X-ray event
          X_ray[p] = 1
      #else:
        #print 'Bad Event!'
        pulse_height[p] = ((f[j+2]&0x3F)<<8) | ((f[j+3]&0xFF))
        j += 4 #move pointer to next event
        p += 1
    else:
        i += 4 + mini_length
    i += 4 + mini_length #move pointer to next minipacket
  return event_time, test, reset, X_ray, pulse_height

#def xray_event_extract():

ph, stat, etime = [], [], []
#read in packets from binary files
LRS, HSK, n, x = telemetry(read_into_buffer(files))
etime = array(etime)
ph = array(ph)
#q = where(array(stat) == 0) # pick X-ray events
q = where((array(stat) == 0)&(ph > 0)) # pick X-ray events
ph = ph[q[0]]*1.0
voltage = 2.5*(ph/(2**14 - 1))
hist, bins = histogram(voltage, bins=1024)

plt.ion()
plt.figure(1)
plt.clf()
plt.plot(bins[1:1025], hist, '-b')
plt.show()

plt.figure(2)
plt.clf()
plt.plot(etime[q])
plt.show()

#list of true X-ray events
# time, tp, rp, X, ph = arange(0, 100, 1) , zeros(100), zeros(100), full(100, 1), zeros(100)
# x, t, r = where(X == 1), where(tp == 1), where(rp == 1)
# #create header for FITS file
# prihdr = pyfits.Header()
# prihdr['MISSION'] = 'HaloSat'
# #date and time of file creation
# prihdr['SYSTIME'] = (str(datetime.now()), 'FITS file creation time')
# #prihdr['MINTIME'] = (time[x].min(), 'Earliest X-ray event time')
# #prihdr['MAXTIME'] = (time[x].max(), 'Latest X-ray event time')
# prihdr['NUMXRAY'] = (len(X[x]), 'Number of X-ray events')
# prihdr['NUMRESET'] = (len(rp[r]), 'Number of reset pulses')
# prihdr['NUMTEST'] = (len(tp[t]), 'Number of test pulses')
# #create primary HDU
# prihdu = pyfits.PrimaryHDU(header=prihdr)
# #create columns for FITS file
# event_time = pyfits.Column(name='Event Time', format='J', array=time)
# test = pyfits.Column(name='Test Event', format='B', array=tp)
# reset = pyfits.Column(name='Reset Pulse', format='B', array=rp)
# X_ray = pyfits.Column(name='X-ray', format='B', array=X)
# pulse_height = pyfits.Column(name='Pulse Height', format='I', array=ph)
# cols = pyfits.ColDefs([event_time, test, reset, X_ray, pulse_height])
# #create binary table HDU
# tbhdu = pyfits.BinTableHDU.from_columns(cols)
# #write header and table to FITS file
# thdulist = pyfits.HDUList([prihdu, tbhdu])
# thdulist.writeto(fits_name)
