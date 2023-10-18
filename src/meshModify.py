import re

def modifyMeshFile( filename ):

    elementStart = False
    val801Starts = False
    val801Done = False
    val812Done = False

    newLine94 = ['94', '3', '2', '15', '15', '7', '5', '9', '8']
    newLine99 = ['98', '3', '2', '15', '15', '14', '8', '15', '12']
    newLine449 = ['447', '3', '2', '15', '15', '342', '344', '346', '345']
    newLine460 = ['457', '3', '2', '15', '15', '346', '354', '353', '355']
    newLine801 = ['797', '3', '2', '15', '15', '667', '668', '671', '670']
    newLine812 = ['807', '3', '2', '15', '15', '671', '679', '678', '680']
    newLine1162 = ['1156', '3', '2', '15', '15', '970', '969', '972', '971']
    newLine1167 = ['1160', '3', '2', '15', '15', '975', '971', '976', '974']

    newLines = []

    with open(filename) as filehandle:

        allLines = filehandle.readlines()
        for lineval in allLines:

            newLines.append( lineval.split() )
            # print(lineval)
            if lineval == "$Elements\n":
                print("found")
                break

        newLines.append(["1167"])

        newsize = len(newLines)
        print(newsize)

    # with open( filename ) as filehandle:

        
        # print(allLines)
        
        for idx, lineval in enumerate(allLines[newsize:]):
            
            if lineval == "$EndElements\n":
                break

            allvals = lineval.split()
            
            if int(allvals[0]) < 94:
                newLines.append( allvals )
                continue

            if int(allvals[0]) == 94:
                newLines.append( newLine94 )
                continue

            if int(allvals[0]) == 95:
                continue

            if int(allvals[0]) == 99:
                newLines.append( newLine99 )
                continue

            if int(allvals[0]) == 103:
                continue

            if int(allvals[0]) == 449:
                newLines.append( newLine449 )
                continue

            if int(allvals[0]) == 450:
                continue

            if int(allvals[0]) == 460:
                newLines.append( newLine460 )
                continue

            if int(allvals[0]) == 462:
                continue

            if int(allvals[0]) == 801:
                newLines.append( newLine801 )
                continue

            if int(allvals[0]) == 802:
                continue

            if int(allvals[0]) == 812:
                newLines.append( newLine812 )
                continue

            if int(allvals[0]) == 814:
                continue

            if int(allvals[0]) == 1162:
                newLines.append( newLine1162 )
                continue

            if int(allvals[0]) == 1163:
                continue

            if int(allvals[0]) == 1167:
                newLines.append( newLine1167 )
                continue

            if int(allvals[0]) == 1170:
                continue

            if int(allvals[0]) > 95 and int(allvals[0]) < 103:
                newLines.append( allvals )
                newLines[-1][0] = str( int( allvals[0] ) - 1 )
                continue

            if int(allvals[0]) > 103 and int(allvals[0]) < 450:
                newLines.append( allvals )
                newLines[-1][0] = str( int( allvals[0] ) - 2 )
                continue

            if int(allvals[0]) > 450 and int(allvals[0]) < 462:
                newLines.append( allvals )
                newLines[-1][0] = str( int( allvals[0] ) - 3 )
                continue

            if int(allvals[0]) > 462 and int(allvals[0]) < 802:
                newLines.append( allvals )
                newLines[-1][0] = str( int( allvals[0] ) - 4 )
                continue

            if int(allvals[0]) > 802 and int(allvals[0]) < 814:
                newLines.append( allvals )
                newLines[-1][0] = str( int( allvals[0] ) - 5 )
                continue

            if int(allvals[0]) > 814 and int(allvals[0]) < 1163:
                newLines.append( allvals )
                newLines[-1][0] = str( int( allvals[0] ) - 6 )
                continue

            if int(allvals[0]) > 1163 and int(allvals[0]) < 1170:
                newLines.append( allvals )
                newLines[-1][0] = str( int( allvals[0] ) - 7 )
                continue

            if int(allvals[0]) > 1170:
                newLines.append( allvals )
                newLines[-1][0] = str( int( allvals[0] ) - 8 )
                continue                
        newLines.append( ["$EndElements"] )

        for idx, lineval in enumerate(newLines):
            newLines[idx] = " ".join(lineval) + '\n'

    newMeshFileName = "/home/gaurav/gmshAutoScripts/build/newMeshModifiedAllCorners.msh"
        
    with open(newMeshFileName, "+w") as filehandle:
        filehandle.writelines( newLines )

    return 0

def modifyMixMeshFile( filename ):

    elementStart = False
    newLines = []

    with open(filename) as filehandle:

        allLines = filehandle.readlines()
        for lineval in allLines:

            newLines.append( lineval.split() )
            # print(lineval)
            if lineval == "$Elements\n":
                print("found")
                break

        newLines.append(["1130"])

        newsize = len(newLines)
        print(newsize)

    # with open( filename ) as filehandle:

        
        # print(allLines)
        
        for idx, lineval in enumerate(allLines[newsize:]):
            
            if lineval == "$EndElements\n":
                break

            allvals = lineval.split()
            
            if int(allvals[0]) < 2:
                newLines.append( allvals )
                continue

            if int(allvals[0]) == 2:
                continue

            if int(allvals[0]) > 2:
                newLines.append( allvals )
                newLines[-1][0] = str( int( allvals[0] ) - 1 )
                continue

        newLines.append( ["$EndElements"] )

        for idx, lineval in enumerate(newLines):
            newLines[idx] = " ".join(lineval) + '\n'

    newMeshFileName = filename
        
    with open(newMeshFileName, "+w") as filehandle:
        filehandle.writelines( newLines )

    return 0

if __name__ == "__main__":

    filename = "/home/gaurav/Finch/src/examples/Mesh/MeshRun/mix_mesh/mesh_lvl4.msh"
    # modifyMeshFile(filename)
    modifyMixMeshFile( filename )
