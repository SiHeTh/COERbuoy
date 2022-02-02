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
    title=header[0];
    var cols=new Array(header.length-stop);
    for (var i=0; i< cols.length; i++){cols[i]={x:[],y:[],name:""}}
    lines.forEach((res,idx)=>{
    var row=res.split(",");
	    if (idx>0)
	    {
		if (row.length == cols.length+stop)
		{for (i=0; i< cols.length; i++){cols[i].x.push(parseFloat(row[0]));cols[i].y.push(parseFloat(row[i+1]));}
		}
	    }
	    else
	    {
	    for (i=0; i< cols.length; i++){cols[i].name=row[i+1]};
	    }
	}
    );
    return {cols, header, title};
}

class ChartFromFile
{
    constructor(type, element)
    {
        this.type=type;
        this.element=document.querySelector(element);
    }
    fromCSV(data,stop,tx,ty,url)
    {
    let csv=readCSV(data,stop);

    
    var series1=[];
    var seriesxy=[];
    for (var i=0; i< csv["cols"].length-1; i++)
    {series1.push(csv["header"][i+1]);seriesxy.push(csv["cols"][i]);}
    //this.chart.destroy();
    this.newChart(series1,seriesxy,tx,ty,csv["title"],url);
    //for (i=0; i< cols.length-1; i++)
    //{this.chart.appendSeries({name:header[i+1],data: cols[i]})}
    //this.chart.resetSeries(true,true);
    return csv["cols"];
    }

    newChart(series1,seriesxy,tx,t,ty,url1)
    {
        if (this.type=='line')
        {
        
        // Define Data
	var data = seriesxy;//[
// 	 {x: x1Values, y: y1Values, mode:"lines"},
// 	 {x: x2Values, y: y2Values, mode:"lines"},
// 	 {x: x3Values, y: y3Values, mode:"lines"}
//	];

	// Define Layout
	var layout = {
	  autosize:"true",
	  xaxis: {title: tx},
	  //yaxis: {title: ty},
	  title:t,
	  height:400,
	  legend:{orientation:"h", y:-0.2, xanchor:'centre'},
	  margin: {
	    l: 0,
	    r: 0,
	    b: 0,
	    t: 30,
	    pad: 1
	  }
	};
		
    console.log(url1)
    let url2=url1;
	 var buttonadd={modeBarButtonsToAdd: [
	    {
	      name: 'Download csv',
	      icon: Plotly.Icons.disk,
	      direction: 'up',
	      click: function(e){window.location.href = url2;}
	      		   
	    }],modeBarButtonsToRemove: ['lasso2d','zoomIn2d','zoomOut2d','resetScale2d']};
    
	// Display using Plotly
	this.chart = Plotly.newPlot(this.element, data, layout,  buttonadd);

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
        this.chart = new ApexCharts(this.element, options);
        this.chart.animations=false;
        this.chart.render();
        }
    
      
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
