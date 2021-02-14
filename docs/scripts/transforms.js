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
    out[ii] = 2 * y0 / (1 + my ** 2) / Math.sqrt(r[ii] ** 2 - y0 ** 2 / (1 + my ** 2)) // forward transform
    if (r[ii]<=(y0 / Math.sqrt(1 + my ** 2))) {out[ii] = 100}
  }
  return out
}
//-------------------------------------------------------------------------//


var inferno = [[0.0952, 0.0459, 0.2380], [0.1502, 0.0453, 0.3185],
  [0.2109, 0.0371, 0.3784], [0.2698, 0.0406, 0.4113],
  [0.3261, 0.0568, 0.4270], [0.3811, 0.0771, 0.4328],
  [0.4358, 0.0974, 0.4320], [0.4906, 0.1170, 0.4256],
  [0.5454, 0.1363, 0.4140], [0.6001, 0.1560, 0.3970],
  [0.6543, 0.1772, 0.3747], [0.7071, 0.2010, 0.3474],
  [0.7578, 0.2287, 0.3157], [0.8055, 0.2614, 0.2802],
  [0.8490, 0.2998, 0.2420], [0.8874, 0.3441, 0.2017],
  [0.9201, 0.3938, 0.1597], [0.9466, 0.4483, 0.1160],
  [0.9667, 0.5067, 0.0705], [0.9802, 0.5683, 0.0301],
  [0.9871, 0.6323, 0.0319], [0.9871, 0.6982, 0.0880],
  [0.9802, 0.7655, 0.1647], [0.9669, 0.8333, 0.2571],
  [0.9509, 0.8992, 0.3717], [0.9500, 0.9562, 0.5112]]
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




var r_vec = linspace(0, 1, 280)
var y0 = 0.5
var ya = abel(y0, r_vec)
var ys = []
var my_vec = linspace(0,4,17)
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
var $container = $('#my_dataviz'),
    width_t = 0.98 * $container.width(),
    height_t = $container.height()
var margin = {
    top: 0,
    right: 10,
    bottom: 50,
    left: 40
  },
  width = width_t - margin.left - margin.right,
  height = 470 - margin.top - margin.bottom;

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
  .call(d3.axisBottom(x).ticks(5))
  .attr("class", "axis");
var xAxis2 = svg.append("g")
  .call(d3.axisTop(x).ticks(0))
  .attr("class", "axis");

// Add Y axis
var y = d3.scaleLinear()
  .domain([0, 5.2])
  .range([height, 0]);
svg.append("g")
  .call(d3.axisLeft(y).ticks(5))
  .attr("class", "axis");
svg.append("g")
  .attr("transform", "translate(" + width + ",0)")
  .call(d3.axisRight(y).ticks(0))
  .attr("class", "axis");

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
  .attr('transform', 'translate(-25,' + height / 2 + ')rotate(-90)')
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
      .transition(100)
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
    .transition(100)
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
