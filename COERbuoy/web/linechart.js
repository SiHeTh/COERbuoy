function readCSV2Table(data,doc)
{	let csvData=[];
    let lines=data.split("\n");
    var header=lines[0].split(",");

    var cols=new Array(header.length);
    for (var i=0; i< cols.length; i++){cols[i]=document.createElement("tr");}
    lines.forEach((res,idx)=>{
    var row=res.split(",");
        if (row.length == cols.length)
        {for (i=0; i< cols.length; i++){var n=document.createElement("th"); n.innerHTML=row[i]; cols[i].append(n)}
        }
    }
    );
    return cols;
}
    
function readCSV(data,stop)
{	let csvData=[];
    let lines=data.split("\n");
    var header=lines[0].split(",");

    var cols=new Array(header.length-stop);
    for (var i=0; i< cols.length; i++){cols[i]=new Array()}
    lines.forEach((res,idx)=>{
    var row=res.split(",");
    if (idx>0)
    {
        if (row.length == cols.length+stop)
        {for (i=0; i< cols.length; i++){cols[i].push({x:parseFloat(row[0]),y:parseFloat(row[i+1])})}
        }}
    });
    return {cols, header};
}

class ChartFromFile
{
    constructor(type, element)
    {
        this.type=type;
        this.element=document.querySelector(element);
    }
    fromCSV(data,stop,tx,ty)
    {
    let csv=readCSV(data,stop);

    
    var series1=[];
    for (var i=0; i< csv["cols"].length-1; i++)
    {series1.push({name:csv["header"][i+1],data: csv["cols"][i]})}
    this.chart.destroy();
    this.newChart(series1,tx,ty);
    //for (i=0; i< cols.length-1; i++)
    //{this.chart.appendSeries({name:header[i+1],data: cols[i]})}
    //this.chart.resetSeries(true,true);
    return csv["cols"];
    }

    newChart(series1,tx,ty)
    {
        if (this.type=='line')
        {
            var options = {
             series: series1
              ,
              chart: {
              
            animations: {
              enabled: false},
              //id: 'line-chart-WEC',
              type: 'line',
              height: 350,
              zoom: {
                autoScaleYaxis: true
              }
            },
            annotations: {
              yaxis: [{
                y: 30,
                borderColor: '#999',
                /*label: {
                  show: true,
                  text: 'Support',
                  style: {
                    color: "#fff",
                    background: '#00E396'
                  }
                }*/
              }],
              xaxis: [{
                //x: 0,
                borderColor: '#999',
                //yAxisIndex: 0,
                /*label: {
                  show: true,
                  text: 'Rally',
                  style: {
                    color: "#fff",
                    background: '#775DD0'
                  }
                }*/
              }]
            },
            dataLabels: {
              enabled: false
            },
            markers: {
              size: 0,
              style: 'hollow',
            },
            xaxis: {
              tickAmount: 6,
              labels:{format:'ff',},
              title:{text:tx}
            },
            yaxis: {
              tickAmount: 6,
              title:{text:ty}
            },
            tooltip: {
              x: {
                format: 'i'
              }
            },
            fill: {
              type: 'solid',
              //gradient: {
                //shadeIntensity: 1,
                //opacityFrom: 0.7,
                //opacityTo: 0.9,
                //stops: [0, 100]
              //}
            },
            };
        }
        else
        {
            var options = {
                series: series1,
                chart: {
                    animations: {
                    enabled: false},
                    height: 350,
                    type: 'heatmap',
                    },
                dataLabels: {
                    enabled: true
                    },
                colors: ["#008FFB"],
                //title: {
                  //  text: 'HeatMap Chart (Single color)'
                    //},
                };
        }
    
        this.chart = new ApexCharts(this.element, options);
        this.chart.animations=false;
        this.chart.render();
      
      }
      setXzoom(a,e)
      {
        this.chart.zoomX(a,e);
      }
      /*  var resetCssClasses = function(activeEl) {
        var els = document.querySelectorAll('button')
        Array.prototype.forEach.call(els, function(el) {
          el.classList.remove('active')
        })
      
        activeEl.target.classList.add('active')
      }
      
      
      document
        .querySelector('#three_waves')
        .addEventListener('click', function(e) {
          resetCssClasses(e)
      
          this.chart.zoomX(
            0,//new Date('27 Sep 2012').getTime(),
            4//new Date('27 Feb 2013').getTime()
          )
        })
      
      document.querySelector('#all').addEventListener('click', function(e) {
        resetCssClasses(e)
      
        this.chart.zoomX(
          0,//new Date('23 Jan 2012').getTime(),
          11//new Date('27 Feb 2013').getTime()
        )
      })
      }
      */
}
