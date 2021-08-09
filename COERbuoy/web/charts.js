function handleFileSelect(evt) {
        var file = evt.target.files[0];
        loadFile(file);}
    function loadFile(file)
    {
        var reader = new FileReader();
        reader.readAsBinaryString(file);
        reader.onload = function(e){data=fromCSV(e.target.result)};
    }
    function fromCSV(data)
    {

    let csvData=[];
    let lines=data.split("\n");
    header=lines[0].split(",");

    var cols=new Array(header.length-1);
    for (i=0; i< cols.length; i++){cols[i]=new Array()}
    lines.forEach((res,idx)=>{
    row=res.split(",");
    if (idx>0)
    {
        if (row.length == cols.length+1)
        {for (i=0; i< cols.length; i++){cols[i].push({x:parseFloat(row[0]),y:parseFloat(row[i+1])})}
        }}
    });

    chart1.destroy();
    newChart();
    for (i=0; i< cols.length-1; i++)
    {chart1.appendSeries({name:header[i+1],data: cols[i]})}
    return cols;
    }

function newChart()
    {
        var options = {
          series: [//{
          //name:"def",
          //data: [
          //]
          //}
          ],
          chart: {
          id: 'line-WEC',
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
            label: {
              show: true,
              text: 'Support',
              style: {
                color: "#fff",
                background: '#00E396'
              }
            }
          }],
          xaxis: [{
            x: 0,
            borderColor: '#999',
            yAxisIndex: 0,
            label: {
              show: true,
              text: 'Rally',
              style: {
                color: "#fff",
                background: '#775DD0'
              }
            }
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
          min: 0,
          tickAmount: 6,
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

        chart1 = new ApexCharts(document.querySelector("#chart-timeline"), options);
        chart1.render();
      
      
        var resetCssClasses = function(activeEl) {
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
      
          chart1.zoomX(
            0,//new Date('27 Sep 2012').getTime(),
            4//new Date('27 Feb 2013').getTime()
          )
        })
      
      document.querySelector('#all').addEventListener('click', function(e) {
        resetCssClasses(e)
      
        chart1.zoomX(
          0,//new Date('23 Jan 2012').getTime(),
          11//new Date('27 Feb 2013').getTime()
        )
      })
      }
