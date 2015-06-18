# This is a script that, upon running, moves the current version of GenSIP to the
# "Previous Versions" directory and makes a new directory with an updated version number. 
# then it iterates through all files and subdirectory files in the new GenSIP package
# and updates the version number therein. 

import os
import sys

cwd = os.getcwd()

directories = os.listdir(cwd)
# Make sure the current version of GenSIP is in the current directory, and assign
# the variable "oldVersion" to the current version's name. 
for f in directories:
    if f.startswith("GenSIP_v"):
        oldVersion = f
        break
else:
    print "There is no GenSIP_vX in current directory."
# Find the oldVersion's version number and create a new version name with the 
# updated number. 
ind = oldVersion.find('v')+1
oldVnum = int(oldVersion[ind:])
newVnum = oldVnum + 1
newName = oldVersion[:ind] + str(newVnum)

# Create a new directory with the updated version number.
if not(os.path.exists(newName)):
    os.mkdir(newName)
    
# Create a list of all files and directories in the old directory 
# Then remove all the irrelevant files: compiled python files and ".DS_Store"
oldFiles = os.listdir(oldVersion)

for filename in oldFiles:
    if filename == ".DS_Store":
        oldFiles.remove(filename)
    # Remove compiled python files (.pyc)
    elif filename.endswith(".pyc"):
        oldFiles.remove(filename)
        
# make a copy of oldFiles as a list called "newFiles."
newFiles = list(oldFiles)

# Create proper paths to old and new files
for filename in oldFiles:
    newFiles[oldFiles.index(filename)] = newName + "/" + filename
    oldFiles[oldFiles.index(filename)] = oldVersion + "/" + filename

# Update the directories
def updateDirectory(oldFiles, newFiles, oldName, newName): 
    # For every file in the old GenSIP folder, create a new file, but replace
    # any instance of the old version name with a new version name.
    for filename in oldFiles:
        newfilename = newFiles[oldFiles.index(filename)]
        # Check if the path is a directory or not. 
        if not(os.path.isdir(filename)):
            if filename.endswith(".py"):
                old = open(filename,'r')
                new = open(newfilename, 'w')
                for line in old:
                    newline = line.replace(oldName, newName)
                    new.write(newline)
                old.close()
                new.close()
            else:
                old = open(filename,'r')
                new = open(newfilename, 'w')
                oldLines = old.readlines()
                for line in oldLines:
                    new.write(line.replace(oldName,newName))
                old.close()
                new.close()
        # if the path is a directory, make a new directory in the new version, and
        # create  a list of contents of the old directory as before, then recursively 
        # apply updateDirectory to the subdirectories. 
        else:
            print "Updating contents of directory %s to %s" % (filename,newfilename)
            os.mkdir(newFiles[oldFiles.index(filename)])
            oldSub = os.listdir(filename)
            if ".DS_Store" in oldSub:
                oldSub.remove(".DS_Store")
            newSub = list(oldSub)
            
            for name in oldSub:
                new = filename + "/" + name
                oldSub[oldSub.index(name)] = new
            for name in newSub:
                new = newfilename + "/" + name
                newSub[newSub.index(name)] = new
            print oldSub
            print newSub
            updateDirectory(oldSub, newSub, oldName, newName)
    
# Run updateDirectory on the old and new files. 
updateDirectory(oldFiles, newFiles, oldVersion, newName)
os.rename(oldVersion, "Previous Versions/"+oldVersion)
print "Finished updating %s to %s" % (oldVersion, newName)
            