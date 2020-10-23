//==========================================================================//
// LINEAR PROJECTION CODE --------------------------------------------------//
var pi = 3.14159265359
var eps = 1e-14 // tolerance

var Kb = function (m, x0, r) {
  return math.re(math.log(math.add(r,
    math.sqrt(math.Complex(r ** 2 - x0 ** 2 / (1 + m ** 2), 0)))))
}
var Kc = function (m, x0, r1, r2, r3) {
  return 1 / (r2 - r1) * Kb(m, x0, r3 + eps)
    // the + eps allows for finite value of kernel when r3 = x0
}
var linear = function (re, y0, my) {
  Nr = re.length - 1;

  var K = Array(re.length)
  for (ii = 0; ii < re.length; ii++) {
    K[ii] = Array(y0.length)
    for (jj = 0; jj < y0.length; jj++) {
      if (ii==0) {
        rjd  = re[ii]; // r_j
        rj = re[ii + 1]; // r_{j+1}
        K[ii][jj] = 2 * y0[jj] * (
          Kc(my[jj],y0[jj],rj,rjd,rj) -
          Kc(my[jj],y0[jj],rj,rjd,rjd)) // decline, first element

      } else if(ii==(re.length - 1)) {
        rj  = re[ii - 1]; // r_j
        rju = re[ii]; // r_{j+1}
        K[ii][jj] = 2 * y0[jj] * (
          Kc(my[jj],y0[jj],rj,rju,rju) -
          Kc(my[jj],y0[jj],rj,rju,rj))

      } else {
        rjd = re[ii - 1]; // r_{j-1}
        rj  = re[ii]; // r_j
        rju = re[ii + 1]; // r_{j+1}
        K[ii][jj] = 2 * y0[jj] * (  // real(.) removes values outside integral bounds
         Kc(my[jj],y0[jj],rjd,rj,rj) -
         Kc(my[jj],y0[jj],rjd,rj,rjd) +  // integral over rise
         Kc(my[jj],y0[jj],rju,rj,rju) -
         Kc(my[jj],y0[jj],rju,rj,rj))
      }

      if (Math.abs(K[ii][jj])<(1e3 * eps)) { K[ii][jj] = 0 } // remove numerical noise

    }
  }
  K = math.transpose(K)
  return K
}

//-------------------------------------------------------------------------//


var r_vec = linspace(0,1,125)

var ya = abel(y0, r_vec)
var ys = []
var my_vec = linspace(0,3,11)
for (jj in my_vec) {
  ys[my_vec[jj].toPrecision(4)] = sipkens(my_vec[jj], y0, r_vec)
}


var re_vec = linspace(0, 1, 50)
var yl = linspace(-1.5, 1.5, 80)
var ml = function (yl, zc, yc) {
  for (ii=0; ii<yl.length; ii++) {
    out[ii] = math.subtract(yl[ii], yc) / zc
  }
  return out
}

var normpdf = function (r, mu, sig) {
  out = Array(r.length)
  for (ii=0; ii<r.length; ii++) {
    out[ii] = 1 / (sig * Math.sqrt(2 * pi)) * Math.exp(-0.5 * ((r[ii] - mu) / sig) ** 2)
    }
  return out
}
var bet1 = normpdf(re_vec, 0, 0.25)
var bet2 = normpdf(re_vec, 0, 0.15)
var bet = Array(re_vec.length)
for (ii=0; ii<re_vec.length; ii++) {
  bet[ii] = 2.2 * bet1[ii] - 1.2 * bet2[ii]
}

var yc = 0
var zc_vec = [10, 6, 3, 2, 1.5, 1.2, 1, 0.8]
var xl = Array(zc_vec.length)

// transfer to data structure
var data2 = [];
for (ii in yl) {
  t0 = {x: yl[ii]}
  for (zz in zc_vec) { // loop through my slopes
    t0['xl' + zz] = 0
  }
  data2.push( t0 )
}

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
var svg2 = d3.select("#my_dataviz2")
  .append("svg")
  .attr("width", width + margin.left + margin.right)
  .attr("height", height + margin.top + margin.bottom)
  .append("g")
  .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

// Add X axis
var x2 = d3.scaleLinear()
  .domain([-1.5, 1.5])
  .range([0, width]);
var xAxis3 = svg2.append("g")
  .attr("transform", "translate(0," + height + ")")
  .call(d3.axisBottom(x2).ticks(5));
var xAxis4 = svg2.append("g")
  .call(d3.axisTop(x2).ticks(0));

// Add Y axis
var y2 = d3.scaleLinear()
  .domain([-7, 16])
  .range([height, 0]);
svg2.append("g")
  .call(d3.axisLeft(y2).ticks(5));
svg2.append("g")
  .attr("transform", "translate(" + width + ",0)")
  .call(d3.axisRight(y2).ticks(0))

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
  .attr('transform', 'translate(-35,' + height / 2 + ')rotate(-90)')
  .text("Projection")

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


// add controls
var updateData2 = function (val) {
  yc = val
  for (zz=0; zz<zc_vec.length; zz++) {
    Kl = linear(re_vec, yl, ml(yl, zc_vec[zz], yc))
    xl[zz] = math.multiply(Kl, bet)
  }

  var data2 = [];
  for (ii in yl) {
    t0 = {x: yl[ii]}
    for (zz in zc_vec) { // loop through my slopes
      t0['xl' + zz] = xl[zz][ii]
    }
    data2.push( t0 )
  }

  for (jj in zc_vec) {
    svg2.select("#p0" + jj)
      .datum(data2)
      .transition(1000)
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

var ycvals = linspace(0, 3, 13)
for (ii in ycvals) { ycvals[ii] = ycvals[ii].toFixed(3) }
function displayycval(val) { // update displayed value
  document.getElementById('ycval').value = ycvals[val - 1];
}
displayycval(document.getElementById('ycSlider').value)
updateData2(ycvals[document.getElementById('ycSlider').value - 1], re_vec) // initial update
d3.select("#ycSlider").on("change", function(d) { // udpate data and plot
  val = this.value
  updateData2(ycvals[val - 1], re_vec)
})
//------------------------------------------------------------------------//
//END PLOT ---------------------------------------------------------------//
