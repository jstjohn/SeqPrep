#!/usr/bin/env python
from glob import iglob as glob
from operator import itemgetter
flist = glob("trimmed_*.txt")
sensdict = dict()
specdict = dict()
xpoints = []
ypoints = []

for fn in flist:
    f = open(fn,'r')
    parts = fn.split("_")
    res_parts = parts[1:-1]
    res_parts.append(parts[-1][:-4]) #remove .txt
    res_parts2 = [ x[0]+":"+x[1:] for x in res_parts]
    rname = '_'.join(res_parts2)
    for line in f:
        if line.startswith("Adapter Trimming Specificity:"):
            spec = float(line.split()[-1])
            xpoints.append(spec)
            specdict[rname] = spec
        elif line.startswith("Adapter Trimming Sensitivity:"):
            sens = float(line.split()[-1])
            ypoints.append(sens)
            sensdict[rname] = sens
    f.close()


outf =  open('spec_sorted_results.txt','w')
colnames = []
rowdata = []
outf.write("#Fname\tSpecificity\tSensitivity\n")
total = len(specdict.items())
num = 1
for (fname,spec) in sorted(specdict.items(),key=itemgetter(1),reverse=True):
    colnames.append("data.addColumn('number','%s');"%(fname))
    sens = sensdict[fname]
    rowarr = [str(sens)]
    for x in range(1,total+1):
        if x == num:
            rowarr.append(str(spec))
        else:
            rowarr.append("null")
    rowdata.append("data.addRow([%s]);" %(", ".join(rowarr)))
    outf.write("%s\t%f\t%f\n"%(fname,spec,sens))
    num += 1
outf.close()

outf = open("result.html","w")


goog_head = '''
<html>
  <head>
    <script type="text/javascript" src="https://www.google.com/jsapi"></script>
    <script type="text/javascript">
      google.load("visualization", "1", {packages:["corechart"]});
      google.setOnLoadCallback(drawVisualization);
      function drawVisualization() {
          // Create and populate the data table.
          var data = new google.visualization.DataTable();
          data.addColumn('number', 'Sensitivity');
          %s
          %s

          // Create and draw the visualization.
          var chart = new google.visualization.ScatterChart(
          document.getElementById('visualization'));
          chart.draw(data, {title: 'Adapter Removal With Different Settings',
              width: 800, height: 800,
              vAxis: {title: "Adapter Specificity",
                  //minValue: 0.8,
                  //maxValue: 1,
                  titleTextStyle: {color: "green"}},
              hAxis: {title: "Adapter Sensitivity",
                  //minValue: 0.8,
                  //maxValue: 1,
                  titleTextStyle: {color: "green"}},          pointSize: 4,       
              legend: 'none'}
           );
      }
    </script>
  </head>

  <body>
    <div id="visualization"></div>
  </body>
</html>
''' % ("\n          ".join(colnames),"\n          ".join(rowdata))

outf.write(goog_head)
outf.close()

