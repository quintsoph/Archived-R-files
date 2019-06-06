####################################################################
import glob
import sys
import statistics
import os
####################################################################


def main():

    #check for arguments by calling parseArgs
    args = parseArgs(sys.argv)
    #if there are no arguments then close program
    if args == ():
        return

    globPattern=args+'/'+'*'
    inputFiles=glob.glob(globPattern)
    for file in inputFiles:
        outFile = open(file[:-4]+"_stripped.txt", 'w')
        outFile.close()
        outFile = open(file[:-4]+"_stripped.txt", 'w')

        with open(file, 'r') as openFile:
            output = []
            for line in openFile:
                line = line.rstrip().split('\t')
                output.append((line[0],line[1]))
        for i in output:
            outFile.write(str(i[0] + '\t' + i[1] + '\n'))
        outFile.close()

    return


'''
Check arguments passed by user and returns arguments.
'''
def parseArgs(args):
    #check for commandline arguments, change compare value to be 1 more than
    #amount you are expecting, if it is less than expected then
    #print out the expected input
    if len(args) != 2:
        print("Usage: THING1 ")
        return ()

    #otherwise set each var to the respected argument
    arg1 = args[1]


    #return all new vars
    return arg1
####################################################################
#call main function to start with
if __name__ == "__main__":
    main()
