from random import randint
from subprocess import check_output
import os
import datetime

class TZ(datetime.tzinfo):
  def utcoffset(self, dt): return datetime.timedelta(minutes=0)
  def dst(self, dt): return None

replacementContainerSymbol = '@'
gmatPath = 'GMAT'


def run(cmd):
  return check_output(cmd, shell = True)


def runScript(templateFileName, params, generatedScriptFileName = None, removeGeneratedScript = True):
  generatedScriptFileName = generateScript(templateFileName, params, generatedScriptFileName)

  command = "{0} --minimize --run --exit \"{1}\"".format(gmatPath, generatedScriptFileName)

  output = run(command)

  if removeGeneratedScript:
    os.remove(generatedScriptFileName)

  return output

def generateScript(templateFileName, params, generatedScriptFileName = None):
  file = open(templateFileName, 'r')
  fileContents = file.read()
  file.close()

  for (key, value) in params.items():
    fileContents = fileContents.replace(replacementContainerSymbol + key + replacementContainerSymbol, str(value))

  pos = fileContents.find('@_for')
  while pos >= 0:
    inPos = fileContents.find('in', pos)
    tagEndPos = fileContents.find('@', inPos)
    endPos = fileContents.find('@_endfor', inPos)
    endTagEndPos = fileContents.find('@', endPos+1)

    varName = fileContents[pos+6:inPos-1]
    arrName = fileContents[inPos+3:tagEndPos]
    loopTemplate = fileContents[tagEndPos+1:endPos]

    loopFinalText = ''
    for arr in params[arrName]:
      stepText = loopTemplate
      for (key, value) in arr.items():
        stepText = stepText.replace(replacementContainerSymbol + varName + '.' + key + replacementContainerSymbol, str(value))

      loopFinalText += stepText

    endTag = fileContents[endPos:endTagEndPos].split(' ')

    if len(endTag) == 2:
      loopFinalText = loopFinalText[0:-int(endTag[1].strip())]

    fileContents = fileContents[0:pos] + loopFinalText + fileContents[endTagEndPos + 1:]

    pos = fileContents.find('@_for')
    

  if generatedScriptFileName is None:
    generatedScriptFileName = 'generated_%i.script' % randint(0, 999999)

  file = open(generatedScriptFileName, 'w')
  file.write(fileContents)
  file.close()

  return generatedScriptFileName

def getReportData(fileName):
  results = []
  file = open(fileName, 'r')
  lines = file.read().strip().split('\n')[1:]
  file.close()

  for line in lines:
    results.append(list(map(float, line.split(' '))))

  return results

def setGmatPath(newPath):
  global gmatPath
  gmatPath = newPath

def getUTCModJulian(dateTime):
  return 21545 + (dateTime - datetime.datetime(2000, 1, 1, 12, 0, 0, 0, tzinfo=TZ())).total_seconds() / 86400