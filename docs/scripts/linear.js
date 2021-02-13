//==========================================================================//
// LINEAR PROJECTION CODE --------------------------------------------------//
var pi = 3.14159265359
var eps = 1e-14 // tolerance

var Kb = function(m, x0, r) {
  return math.re(math.log(math.add(r,
    math.sqrt(math.Complex(r ** 2 - x0 ** 2 / (1 + m ** 2), 0)))))
}
var Kc = function(m, x0, r1, r2, r3) {
  return 1 / (r2 - r1) * Kb(m, x0, r3 + eps)
  // the + eps allows for finite value of kernel when r3 = x0
}
var linear = function(re, y0, my) {
  var K = Array(re.length)
  for (ii = 0; ii < re.length; ii++) {
    K[ii] = Array(y0.length)
    for (jj = 0; jj < y0.length; jj++) {
      if (ii == 0) {
        rjd = re[ii]; // r_j
        rj = re[ii + 1]; // r_{j+1}
        K[ii][jj] = 2 * y0[jj] / (1 + my[jj] ** 2) * (
          Kc(my[jj], y0[jj], rj, rjd, rj) -
          Kc(my[jj], y0[jj], rj, rjd, rjd)) // decline, first element

      } else if (ii == (re.length - 1)) {
        rj = re[ii - 1]; // r_j
        rju = re[ii]; // r_{j+1}
        K[ii][jj] = 2 * y0[jj] / (1 + my[jj] ** 2) * (
          Kc(my[jj], y0[jj], rj, rju, rju) -
          Kc(my[jj], y0[jj], rj, rju, rj))

      } else {
        rjd = re[ii - 1]; // r_{j-1}
        rj = re[ii]; // r_j
        rju = re[ii + 1]; // r_{j+1}
        K[ii][jj] = 2 * y0[jj] / (1 + my[jj] ** 2) * ( // real(.) removes values outside integral bounds
          Kc(my[jj], y0[jj], rjd, rj, rj) -
          Kc(my[jj], y0[jj], rjd, rj, rjd) + // integral over rise
          Kc(my[jj], y0[jj], rju, rj, rju) -
          Kc(my[jj], y0[jj], rju, rj, rj))
      }

      if (Math.abs(K[ii][jj]) < (1e3 * eps)) {
        K[ii][jj] = 0
      } // remove numerical noise

    }
  }
  K = math.transpose(K)
  return K
}

//-------------------------------------------------------------------------//


var re_vec = linspace(0, 1, 100)
var yl = linspace(-1.25, 1.25, 100)
var ml = function(yl, zc, yc) {
  for (ii = 0; ii < yl.length; ii++) {
    out[ii] = math.subtract(yl[ii], yc) / zc
  }
  return out
}

var normpdf = function(r, mu, sig) {
  out = Array(r.length)
  for (ii = 0; ii < r.length; ii++) {
    out[ii] = 1 / (sig * Math.sqrt(2 * pi)) * Math.exp(-0.5 * ((r[ii] - mu) / sig) ** 2)
  }
  return out
}
var bet1 = normpdf(re_vec, 0, 0.25)
var bet2 = normpdf(re_vec, 0, 0.15)
var bet = Array(re_vec.length)
for (ii = 0; ii < re_vec.length; ii++) {
  bet[ii] = 2.2 * bet1[ii] - 1.2 * bet2[ii]
}

var yc = 0
var zc_vec = [20, 4, 2.5, 2, 1.6, 1.35, 1.2, 1.1, 1, 0.9, 0.8]  // curently unused

var zc = 1
var yc_vec = linspace(3, 0, 11)
var xl = Array(yc_vec.length)

// transfer to data structure
var data2 = []
for (ii in yl) {
  t0 = {
    x: yl[ii]
  }
  for (yy in zc_vec) { // loop through my slopes
    t0['xl' + yy] = 0
  }
  data2.push(t0)
}

data3 = []
for (ii=(re_vec.length-1); ii>=0; ii--) {
  data3.push( { x: -re_vec[ii], y: bet[ii] } )
}
for (ii=0; ii<re_vec.length; ii++) {
  data3.push( { x: re_vec[ii], y: bet[ii] } )
}

//-------------------------------------------------------------------------//
// GENERATE PLOTS
// set the dimensions and margins of the graph
// use height / width / padding from previous plot

// UPPER PANEL ------------------------------------------------------------//
// append the svg object to the body of the page
var $container = $('#my_beta'),
    width_a = 0.98 * $container.width(),
    height_a = $container.height()
var margin = {
    top: 0,
    right: 5,
    bottom: 6,
    left: 45
  },
  width = width_a - margin.left - margin.right,
  height3 = 120 - margin.top - margin.bottom;

var svg3 = d3.select("#my_beta")
  .append("svg")
  .attr("width", width + margin.left + margin.right)
  .attr("height", height3 + margin.top + margin.bottom)
  .append("g")
  .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

// Add X axis
var x2 = d3.scaleLinear()
  .domain([-1.25, 1.25])
  .range([0, width]);
var xAxis2a = svg3.append("g")
  .attr("transform", "translate(0," + height3 + ")")
  .call(d3.axisBottom(x2).ticks(5))
  .attr("class", "axis");
var xAxis2b = svg3.append("g")
  .call(d3.axisTop(x2).ticks(0))
  .attr("class", "axis");

// Add Y axis
var y3 = d3.scaleLinear()
  .domain([0, 1.65])
  .range([height3, 0]);
svg3.append("g")
  .call(d3.axisLeft(y3).ticks(4))
  .attr("class", "axis");
svg3.append("g")
  .attr("transform", "translate(" + width + ",0)")
  .call(d3.axisRight(y3).ticks(0))
  .attr("class", "axis");
// Y axis label:
svg3.append("text")
  .attr("text-anchor", "middle")
  .attr('transform', 'translate(-32,' + height3 / 2 + ')rotate(-90)')
  .html("δ = 1 - n/n0")

// add phantom to panel
svg3.append("path")
  .datum(data3)
  .attr("id", "betpath")
  .attr("fill", "none")
  .attr("stroke", "black")
  .attr("stroke-width", 2)
  .attr("d", d3.line()
    .x(function(d) {
      return x2(d.x)
    })
    .y(function(d) {
      return y3(d.y)
    })
  )





// MAIN PLOT ---------------------------------------------------------------//
// append the svg object to the body of the page
var svg2 = d3.select("#my_dataviz2")
  .append("svg")
  .attr("width", width_a + margin.left + margin.right)
  .attr("height", height + margin.top + margin.bottom)
  .append("g")
  .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

// Add X axis
// inherit x scale from above
var xAxis2a = svg2.append("g")
  .attr("transform", "translate(0," + height + ")")
  .call(d3.axisBottom(x2).ticks(5))
  .attr("class", "axis");
var xAxis2b = svg2.append("g")
  .call(d3.axisTop(x2).ticks(0))
  .attr("class", "axis");

// Add Y axis
var y2 = d3.scaleLinear()
  .domain([-5, 5])
  .range([height, 0]);
svg2.append("g")
  .call(d3.axisLeft(y2).ticks(5))
  .attr("class", "axis");
svg2.append("g")
  .attr("transform", "translate(" + width + ",0)")
  .call(d3.axisRight(y2).ticks(0))
  .attr("class", "axis");

//-- Add axis labels --//
// Add X axis label:
svg2.append("text")
  .attr("text-anchor", "middle")
  .attr('x', width / 2)
  .attr('y', height + 35)
  .text("Radial position, y0 [a.u.]");

// Y axis label:
svg2.append("text")
  .attr("text-anchor", "middle")
  .attr('transform', 'translate(-28,' + height / 2 + ')rotate(-90)')
  .text("Projected deflection at z = 0 [a.u.]")

// generate plot
for (jj in zc_vec) {
  svg2.append("path")
    .datum(data2)
    .attr("id", "p0" + jj)
    .attr("fill", "none")
    .attr("stroke", cInterp(jj, 0, zc_vec.length - 1, inferno))
    .attr("stroke-width", 2)
    .attr("d", d3.line()
      .x(function(d) {
        return x2(d.x)
      })
      .y(function(d) {
        return y2(d['xl' + jj])
      })
      // .defined(((d, i) => d['ys' + jj] != null))
    )
}


// Update plots (for zc_vec), currently unused
var updateData2 = function(valy) {
  yc = valy
  for (zz = 0; zz < zc_vec.length; zz++) {
    Kl = linear(re_vec, yl, ml(yl, zc_vec[zz], yc))
    xl[zz] = math.multiply(Kl, bet)
  }

  var data2 = [];
  for (ii in yl) {
    t0 = {
      x: yl[ii]
    }
    for (zz in zc_vec) { // loop through my slopes
      t0['xl' + zz] = xl[zz][ii]
    }
    data2.push(t0)
  }

  for (jj in zc_vec) {
    svg2.select("#p0" + jj)
      .datum(data2)
      .transition(150)
      .attr("d", d3.line()
        .x(function(d) {
          return x2(d.x)
        })
        .y(function(d) {
          return y2(d['xl' + jj])
        })
        // .defined(((d, i) => d['ys' + jj] != null))
      )
  }
}

// Update plots (for yc_vec)
var updateData3 = function(valz) {
  zc = valz
  for (yy = 0; yy < yc_vec.length; yy++) {
    Kl = linear(re_vec, yl, ml(yl, zc, yc_vec[yy]))
    xl[yy] = math.multiply(Kl, bet)
  }

  var data2 = [];
  for (ii in yl) {
    t0 = {
      x: yl[ii]
    }
    for (yy in yc_vec) { // loop through my slopes
      t0['xl' + yy] = xl[yy][ii]
    }
    data2.push(t0)
  }

  for (jj in yc_vec) {
    svg2.select("#p0" + jj)
      .datum(data2)
      .transition(150)
      .attr("d", d3.line()
        .x(function(d) {
          return x2(d.x)
        })
        .y(function(d) {
          return y2(d['xl' + jj])
        })
        // .defined(((d, i) => d['ys' + jj] != null))
      )
  }
}

updateData2(document.getElementById('zcSlider').value, re_vec) // initial update
d3.select("#zcSlider").on("change", function(d) { // udpate data and plot
  val = this.value
  document.getElementById('zcSlider').disabled = true  // temporarily disable during update
  updateData3(val, re_vec)
  document.getElementById('zcSlider').disabled = false  // re-enable
})
//------------------------------------------------------------------------//
//END PLOT ---------------------------------------------------------------//
