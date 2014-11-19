#!/usr/bin/python

import sys
#import linecache
#import time
#import MySQLdb
import pandas as pd


def load_params(paramsfile, indexer='subid'):
	"""opens params file and creates pandas dataframe
		by default it uses the variable "subid" as the index"""
	paramsframe=pd.read_csv(paramsfile)
	return paramsframe





def load (conf):
    global dict
    dict = {}
    try:
	try:
	    fd = open ( conf, "r" )
	except:
	    fd = open ( "./session.config", "r" )
    except IOError, (errorno,strerror):
	print strerror

    while 1:
	  line = fd.readline ()
	  if not line		      : break
	  if line[0] == "#"	      : continue
	  if len(line.split("=")) < 2 : continue 
	  
	  key	= line.split('=') [0]
	  value = line.split('=') [1]

	  dict [key.strip()]=value.strip()

    fd.close ()

    
#...................db access.................................................

def run (db,query):
    #print "running query.............."+query
    db.query ( query )
    res = db.use_result ()
    rows = []
    while ( res ):
        cols = []
        row = res.fetch_row  ()
        if  len(row) == 0:
            break;
        for i in range(len(row[0])):
            cols.append ( row[0][i] )
        rows.append(cols)
    return rows

def connect(config):
    db_name = get('MYSQLDB',config)
    host = get('MYSQLHOST',config)
    user = get('MYSQLUSER',config)
    passwd = get('MYSQLPASSWD',config)
    #print "connecting....."
    #print "HOST ::::: "+host+" USER ::::: "+user+" DATABASE ::::: "+db_name
    try:
        g_db = MySQLdb.connect( host , user , passwd , db_name );

    except MySQLdb.Error,e:
        print 'Fatal Error ... '
        print 'Error connecting to database ...\n', e.args[0], e.args[1]
        sys.exit (-1)

    db = g_db
  #  print "connected"
#    run( db,"set autocommit=1" )
    return db

# return a single item from the database

def getItem( query,config ):
    item = "None"
    conn = connect(config)
    res = run(conn,query)
    print "res is "+str(res)+" LEN IS "+str(len(res))
    if len(res) > 0:
        item = res[0][0]
    return str(item)
    

#..............................................................................

def get (key,conf):
    global dict
    #if not dict: load (conf)
    load (conf)

    try:
        return dict [key]
    except:
        return None

#..............................................................................

# return floating point for comparing times
    
def getTime(datestr,timestr):
    t = timestr.split('.')[0].strip()
    fulltime = datestr+t
    t_tuple = time.strptime(fulltime, "%Y%m%d%H%M%S")
    t_float = time.mktime(t_tuple)
    return t_float
    
def isTrigger(nums):
    if nums.find("=")>-1:
        return False
    if len(nums.split())<5:
        return False
    else:
        t = float(nums.split()[4].strip())
        if t>3.5:
            return True
    return False
    
def hasTriggers(fname,currLine):
    lineNum = currLine + 7
    triggers = False
    eof = False
    while not eof and not triggers:
        ln = linecache.getline(fname,lineNum)
        if ln == '':
            eof = True
        elif ln.find("=")>-1:
            eof = True
        else:
            triggers = isTrigger(ln)
        lineNum = lineNum + 1
    return triggers

def countlines(fname):
    fd = open(fname,"r")
    eof = False
    lncnt = 0
    while not eof:
        ln = fd.readline()
        if ln != '':
            lncnt = lncnt + 1
        else:
            eof = True
    fd.close()
    return lncnt

    




    
