//==========================================================================//
// TRANSFORMS CODE ---------------------------------------------------------//
// create an array with linear spacing
var linspace = function(a, b, n) {
  if (typeof n === "undefined") n = Math.max(Math.round(b - a) + 1, 1);
  if (n < 2) {
    return n === 1 ? [a] : [];
  }
  var i, ret = Array(n);
  n--;
  for (i = n; i >= 0; i--) {
    ret[i] = (i * b + (n - i) * a) / n;
  }
  return ret;
};

// abel transform
var abel = function (y0, r) {
  out = Array(r.length)
  for (ii in r) {
    out[ii] = 2 * y0 / Math.sqrt(r[ii] ** 2 - y0 ** 2) // forward transform
    if (r[ii]<=y0) {out[ii] = 100}
  }
  return out
}

// new transform from Sipkens et al.
var sipkens = function (my, y0, r) {
  out = Array(r.length)
  for (ii in r) {
    out[ii] = 2 * y0 / Math.sqrt(r[ii] ** 2 - y0 ** 2 / (1 + my ** 2)) // forward transform
    if (r[ii]<=(y0 / Math.sqrt(1 + my ** 2))) {out[ii] = 100}
  }
  return out
}
//-------------------------------------------------------------------------//


var inferno = [[0.0676, 0.0390, 0.1893],
  [0.1285, 0.0473, 0.2896],
  [0.2008, 0.0380, 0.3704],
  [0.2722, 0.0411, 0.4123],
  [0.3399, 0.0618, 0.4292],
  [0.4062, 0.0864, 0.4331],
  [0.4723, 0.1105, 0.4283],
  [0.5386, 0.1339, 0.4157],
  [0.6047, 0.1577, 0.3953],
  [0.6698, 0.1838, 0.3672],
  [0.7328, 0.2143, 0.3321],
  [0.7920, 0.2513, 0.2909],
  [0.8455, 0.2964, 0.2452],
  [0.8918, 0.3500, 0.1965],
  [0.9296, 0.4115, 0.1454],
  [0.9583, 0.4795, 0.0914],
  [0.9774, 0.5526, 0.0380],
  [0.9869, 0.6296, 0.0305],
  [0.9864, 0.7093, 0.0998],
  [0.9759, 0.7909, 0.1973],
  [0.9573, 0.8723, 0.3207]]
function componentToHex(c) {
  var hex = c.toString(16);
  return hex.length == 1 ? "0" + hex : hex;
}
function rgbToHex(r, g, b) {
  return "#" + componentToHex(Math.round(r * 255)) +
  componentToHex(Math.round(g * 255)) +
  componentToHex(Math.round(b * 255));
}
function cInterp(no, min, max, cm) {
  idx1 = (no - min) / (max - min)
  idx2 = idx1 * (cm.length - 2)
  a0 = Math.floor(idx2)
  r = (idx2 - a0) * cm[a0 + 1][0] + (1 - (idx2 - a0)) * cm[a0][0] // interpolated red
  g = (idx2 - a0) * cm[a0 + 1][1] + (1 - (idx2 - a0)) * cm[a0][1]
  b = (idx2 - a0) * cm[a0 + 1][2] + (1 - (idx2 - a0)) * cm[a0][2]
  hex = rgbToHex(r, g, b) // convert to hex code
  return hex
}




var r_vec = linspace(0,1,125)
var y0 = 0.5
var ya = abel(y0, r_vec)
var ys = []
var my_vec = linspace(0,3,11)
for (jj in my_vec) {
  ys[my_vec[jj].toPrecision(4)] = sipkens(my_vec[jj], y0, r_vec)
}

// transfer to data structure
var data = [];
for (ii in r_vec) {
  t0 = {x: r_vec[ii], ya: ya[ii]}
  for (jj in my_vec) { // loop through my slopes
    t0['ys' + jj] = ys[my_vec[jj].toPrecision(4)][ii]
  }
  data.push( t0 )
}
console.log('data = ')
console.log(data)

//-------------------------------------------------------------------------//
// GENERATE PLOTS
// set the dimensions and margins of the graph
var margin = {
    top: 0,
    right: 1.5,
    bottom: 50,
    left: 45
  },
  width = 450 - margin.left - margin.right,
  height = 350 - margin.top - margin.bottom;

// append the svg object to the body of the page
var svg = d3.select("#my_dataviz")
  .append("svg")
  .attr("width", width + margin.left + margin.right)
  .attr("height", height + margin.top + margin.bottom)
  .append("g")
  .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

var svg2 = d3.select("#my_labels")
  .append("svg")
  .attr("width", width + margin.left + margin.right)
  .attr("height", 40)
  .append("g")
  .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

// Add X axis
var x = d3.scaleLinear()
  .domain([0, 1])
  .range([0, width]);
var xAxis = svg.append("g")
  .attr("transform", "translate(0," + height + ")")
  .call(d3.axisBottom(x).ticks(5));
var xAxis2 = svg.append("g")
  .call(d3.axisTop(x).ticks(0));

// Add Y axis
var y = d3.scaleLinear()
  .domain([0, 10.2])
  .range([height, 0]);
svg.append("g")
  .call(d3.axisLeft(y).ticks(5));
svg.append("g")
  .attr("transform", "translate(" + width + ",0)")
  .call(d3.axisRight(y).ticks(0))

//-- Add axis labels --//
// Add X axis label:
svg.append("text")
  .attr("text-anchor", "middle")
  .attr('x', width / 2)
  .attr('y', height + 35)
  .text("Radius, r [a.u.]");

// Y axis label:
svg.append("text")
  .attr("text-anchor", "middle")
  .attr('transform', 'translate(-35,' + height / 2 + ')rotate(-90)')
  .text("Transform")

// generate plot
for (jj in my_vec) {
  svg.append("path")
    .datum(data)
    .attr("id", "sl" + jj)
    .attr("fill", "none")
    .attr("stroke", cInterp(jj, 0, my_vec.length - 1, inferno))
    .attr("stroke-width", 2)
    .attr("d", d3.line()
      .x(function(d) {
        return x(d.x)
      })
      .y(function(d) {
        return y(d['ys' + jj])
      })
      // .defined(((d, i) => d['ys' + jj] != null))
    )
  }

// add Abel curve (dashed)
svg.append("path")
  .datum(data)
  .attr("id", "al")
  .attr("fill", "none")
  .attr("stroke", "white")
  .attr('stroke-dasharray', "4 3")
  .attr("stroke-width", 0.6)
  .attr("d", d3.line()
    .x(function(d) {
      return x(d.x)
    })
    .y(function(d) {
      return y(d.ya)
    })
  )


// add controls
var updateData = function (y0, r_vec) {

  ya = abel(y0, r_vec)
  ys = []
  for (jj in my_vec) {
    ys[my_vec[jj].toPrecision(4)] = sipkens(my_vec[jj], y0, r_vec)
  }

  data = [];
  for (ii in r_vec) {
    t0 = {x: r_vec[ii], ya: ya[ii]}
    for (jj in my_vec) { // loop through my slopes
      t0['ys' + jj] = ys[my_vec[jj].toPrecision(4)][ii]
    }
    data.push( t0 )
  }

  for (jj in my_vec) {
    svg.select("#sl" + jj)
      .datum(data)
      .transition(1000)
      .attr("d", d3.line()
        .x(function(d) {
          return x(d.x)
        })
        .y(function(d) {
          return y(d['ys' + jj])
        })
      )
    }
  // Abel curve (dashed)
  svg.select("#al")
    .datum(data)
    .transition(1000)
    .attr("d", d3.line()
      .x(function(d) {
        return x(d.x)
      })
      .y(function(d) {
        return y(d.ya)
      })
    )
}

var y0vals = linspace(0, 0.9, 19)
for (ii in y0vals) { y0vals[ii] = y0vals[ii].toFixed(3) }
function displayy0val(val) { // update displayed value
  document.getElementById('y0val').value = y0vals[val - 1];
}
displayy0val(document.getElementById('y0Slider').value)
d3.select("#y0Slider").on("change", function(d) { // udpate data and plot
  val = this.value
  updateData(y0vals[val - 1], r_vec)
})
//------------------------------------------------------------------------//
//END PLOT ---------------------------------------------------------------//
