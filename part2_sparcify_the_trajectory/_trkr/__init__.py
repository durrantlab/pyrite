import os
import urllib2
import sys
import urllib

def trkr():
    # Get the user name
    if os.path.exists(".trkr.dat"):
        username = open(".trkr.dat").read().strip()
    else:
        username = raw_input("Please enter your durrantlab.com username: ")
        open(".trkr.dat", 'w').write(username)
    
    # Inform server
    url = 'https://durrantlab.com/usr/trkr.php?username='
    url += urllib.quote_plus(username)
    url += "&msg="
    url += urllib.quote_plus("Using " + os.path.basename(sys.argv[0]))
    response = urllib2.urlopen(url)
    