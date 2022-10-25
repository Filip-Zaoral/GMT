
def write(Log,name):
    
    # Write the message console outputs to the Log File:
    logFile = open(name + '.log','wt')
    logFile.writelines('\n'.join(Log))
    logFile.close()                                                            # The saved file can be located in the working directory.
    return
