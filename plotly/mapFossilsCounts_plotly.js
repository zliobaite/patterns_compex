// var locfile = 'data/resorted_localities_more.csv';
// var xpsfld = 'data/xps/';

// COLORMAPPING
//--------------------------
function enforceBounds(x) {
    if (x < 0) {
        return 0;
    } else if (x > 1){
        return 1;
    } else {
        return x;
    }
}

function interpolateLinearly(x, values) {

    // Split values into four lists
    var x_values = [];
    var r_values = [];
    var g_values = [];
    var b_values = [];
    for (i in values) {
        x_values.push(values[i][0]);
        r_values.push(values[i][1][0]);
        g_values.push(values[i][1][1]);
        b_values.push(values[i][1][2]);
    }

    var i = 1;
    while (x_values[i] < x) {
        i = i+1;
    }
    i = i-1;

    var width = Math.abs(x_values[i] - x_values[i+1]);
    var scaling_factor = (x - x_values[i]) / width;

    // Get the new color values though interpolation
    var r = r_values[i] + scaling_factor * (r_values[i+1] - r_values[i])
    var g = g_values[i] + scaling_factor * (g_values[i+1] - g_values[i])
    var b = b_values[i] + scaling_factor * (b_values[i+1] - b_values[i])

    return "rgb("+(255*enforceBounds(r)).toFixed(0)+","+(255*enforceBounds(g)).toFixed(0)+","+(255*enforceBounds(b)).toFixed(0)+")";

}

// SIMPLE TOOLS
//---------------------------
function unionList(listA, listB) {
    let _union = new Set(listA.concat(listB))
    return [..._union];
}
function dupliList(listA) {
    // return [true].concat(listA, listA, listA, listA, truesList(listA.length-1, 1));
    return listA.concat(listA);
}
function truesList(n, k=2) {
    truesL = []; 
    for (var i = 0; i < k*n; i++) {
        truesL.push(Boolean(1));
    }
    return truesL;
}

function getAssociatedIds(i, n) {
    let bid = -1; 
    if (i >= 0 && i < n){
        bid = i;
    }
    if (i >= n && i < 2*n){
        bid = i-n;
    }
    // console.log(i,n, bid, bid+n);
    if (bid > -1){
        return [bid, bid+n];
    }
    return [];

}

function add(s, v) {
    return s+v
}

function cart2pol(x, y){
    return {rho: Math.sqrt(x*x + y*y),
            phi: Math.atan2(y, x)}
}
    
function pol2cart(rho, phi){
    return {x: rho * Math.cos(phi),
            y: rho * Math.sin(phi)}
}
function scaleVal(v, fact=100, off=0, pw=2){
    // console.log(v);
    if (pw == 1){
        return off-fact*v;
    }
    return off-fact*Math.pow(v, pw);
}

// SELECTORS TOOLS
//---------------------------
function clearOptions(selector) {
    while(selector.firstChild) selector.removeChild(selector.lastChild);
}

function assignOptions(textArray, selector) {
    for (var i = 0; i < textArray.length; i++) {
        var currentOption = document.createElement('option');
        currentOption.text = textArray[i];
        selector.appendChild(currentOption);
    }
}

// DATA TOOLS
//---------------------------
function rtrnLoc(d, v_count, v_focus) {
    return {
        lidnum: +d.LIDNUM, // convert to number
        lat: +d.LAT, // convert to number
        lng: +d.LONG, // convert to number
        maxage: +d.MAX_AGE, // convert to number
        minage: +d.MIN_AGE, // convert to number
        slice: +d.SLICE_ID, // convert to number
        vcount: +d[v_count], // convert to number
        vfocus: +d[v_focus] // convert to number
    };
}
function prepareLocs(data, noiseLoc) {
    var locs = {lidnum: [], lat: [], lng: [], maxage: [], minage: [], slice: [], vcount: [], vfocus: []};
    for (var i = 1; i < data.length; i++) {
        var datum = data[i];
        locs.lidnum.push(datum.lidnum);
        locs.lat.push(datum.lat + noiseLoc*Math.random());
        locs.lng.push(datum.lng + noiseLoc*Math.random());
        locs.maxage.push(datum.maxage);
        locs.minage.push(datum.minage);
        locs.slice.push(datum.slice);
        if ( isNaN(datum.vcount) ){
            locs.vcount.push(-1);
        }
        else {
            locs.vcount.push(datum.vcount);
        }
        if ( isNaN(datum.vfocus) ){
            locs.vfocus.push(-1);
        }
        else {
            locs.vfocus.push(datum.vfocus);
        }
    };
    return locs
}

function prepareLocsColors(dataClass, cparameters) {
    var locsColors = [];
    var v_min = cparameters["v_min"];
    var v_max = cparameters["v_max"];
    if ( v_max < 0 ){
            v_max = Math.max.apply(Math, dataClass);        
    }
    var v_denom = (v_max-v_min);
    if ( v_denom <= 0 ){
        v_denom = 1;
    }   
    for (var i = 0; i < dataClass.length; i++) {
        if (dataClass[i] < v_min){            
            locsColors[i] = "rgb(200,200,200)"; 
        }  else {
            locsColors[i] = interpolateLinearly((dataClass[i]-v_min)/v_denom, ccolormaps[cparameters["cmap"]]);
        }
        // if (dataClass[i] > 0){            
        //     locsColors[i] = "#A1BE56"; //
        // } else {
        //     locsColors[i] = "#781C81";
        // }

    }
    return locsColors;
}

function prepareAges(locs){    
    ages = unionList(locs.maxage, locs.minage).sort(function(a, b){return b-a});
    var agesMap = {};
    for (age of ages) {
        agesMap[age] = {vect: [], count: 0};
        for (var i = 0; i < locs.minage.length; i++) {
            agesMap[age].vect.push(Boolean((locs.minage[i] < age) & (locs.maxage[i] >= age)));
        }
        agesMap[age].count = agesMap[age].vect.reduce(add, 0);        
    }
    return [ages, agesMap]
}
function pairCmp(pairA, pairB){
    if ( pairA[0] == pairB[0] ){
        if ( pairA[1] == pairB[1] ){
            if ( pairA[2] == pairB[2] ){
                if ( pairA[3] == pairB[3] ){
                    return pairB[4] - pairA[4];
                }
                else {
                    return pairA[3] - pairB[3];
                }
            }
            else {
                return pairB[2] - pairA[2];
            }
        }
        else {
            return pairB[1] - pairA[1];
        }
    }
    else {
        return pairA[0] - pairB[0];
    }
}
function idsSorted(vA, vB, vC, vD, vE, vX){
    let ids = [...Array(vA.length).keys()].sort(function(a, b){return pairCmp([vX[a], vA[a], vB[a], vC[a], vD[a], vE[a]], [vX[b], vA[b], vB[b], vC[b], vD[b], vE[b]])});
    var idsMap = {};
    var off = 0;
    var step = 1;
    // var idsK = {}; // to get the sorted list of ids
    for (var i = 0; i < ids.length; i++) {
        if (i > 0  & vX[ids[i-1]] != vX[ids[i]]) {
            off = 0;
            step = -1*step;
        }
        off = off+step;
        idsMap[ids[i]] = off;
        // idsK[i] = vE[ids[i]];
    }
    // console.log(idsK);    
    return idsMap;
}

// PLOTTING TOOLS
//---------------------------
function makeSlider(ages, agesMap){
    var sliderSteps = [];
    sliderSteps.push({
        method: 'update',
        label: "All",
        args: [{'visible': truesList(agesMap[age].vect.length)}]
        });
    for (age of ages) {
        sliderSteps.push({
            method: 'update',
            label: age.toFixed(4),
            args: [{'visible': dupliList(agesMap[age].vect)}]
        });
    }
    return sliderSteps;
}

function makeTraces(locs, ordAges, locsColors, cparameters){
    var frng = 0.012*Math.min(cparameters["xlims"][1]-cparameters["xlims"][0], cparameters["ylims"][1]-cparameters["ylims"][0]);
    idsMapAge = idsSorted(locs.maxage, locs.minage, locs.lat, locs.lng, locs.lidnum, locs.slice);
    var traces = {"mrk": [], "time": []}
    for (var i = 0; i < locs.lidnum.length; i++) {
        let xy_pin = {x: locs.lng[i], y: locs.lat[i]}
        let cc = '#ffffff';
        traces["mrk"].push({lon: [xy_pin.x], 
                            lat: [xy_pin.y],                            
                            text: [""+locs.lidnum[i]+" ("+cparameters["slice_names"][locs.slice[i]]+") "+cparameters["v_focus"]+"="+locs.vfocus[i]],
                            //+"\n["+locs.maxage[i].toFixed(4)+","+locs.minage[i].toFixed(4)+"]"],
                            visible: true,
                            type:'scattergeo',
                            mode:'markers',
                            name: "dot"+i,
                            marker: {
                                // opacity : .66,
                                color: [locsColors[i]],
                                // size: baseMrkSz+Math.sqrt(locs.vcount[i])
                                size: cparameters["baseMrkSz"]+Math.log2(2+locs.vcount[i])
                            },
                            hoverinfo: "text",
                            showlegend: false
                           });
        traces["time"].push({x: [ordAges[locs.maxage[i]]+cparameters["padT"], ordAges[locs.minage[i]]-cparameters["padT"]], 
                             y: [idsMapAge[i], idsMapAge[i]],
                             text: [""+locs.lidnum[i]+" ("+locs.lat[i].toFixed(2)+","+locs.lng[i].toFixed(2)+") "+cparameters["v_focus"]+"="+locs.vfocus[i]],
                             visible: true,
                             type:'scatter',
                             mode:'lines',
                             name: "time"+i,                             
                             line: {
                                 color: locsColors[i],
                                 width: cparameters["baseLineW"]
                             },
                             hoverinfo: "text",
                             showlegend: false,
                             xaxis: 'x2',
                             yaxis: 'y2'
                        });

    };
    return traces;
}

function makePlot(dataLocs, cparameters){    
    var fdims = cparameters["fdims"];

    var baseLineW = cparameters["baseLineW"];
    var baseMrkSz = cparameters["baseMrkSz"];
    var noiseLoc = cparameters["noiseLoc"];
    var padT = cparameters["padT"];

    var v_name = cparameters["v_name"];
    
    var xlims = cparameters["xlims"];
    var ylims = cparameters["ylims"];
    var ttickTxts = cparameters["ttickTxts"];
    var agesTicks = cparameters["agesTicks"];

    // prepare data
    locs = prepareLocs(dataLocs, noiseLoc);
    [ages, agesMap] = prepareAges(locs);

    var ordAges = {}
    for (var i = 0; i < ages.length; i++) {
        xx = (ages[i]-ages[0])/(ages[ages.length-1]-ages[0]);
        ordAges[ages[i]] = xx;
        // tticks.push(xx);
        // ordAges[ages[i]] = i+1;
        // tticks.push(i+1);
        // ttickTxts.push(ages[i].toFixed(4));
    }
    var tticks = [-padT];
    for (var i = 0; i < agesTicks.length; i++) {
        tticks.push((agesTicks[i]-ages[0])/(ages[ages.length-1]-ages[0]));
    }
    tticks.push(1+padT);
    locsColors = prepareLocsColors(locs.vfocus, cparameters);
    
    // prepare plot elements
    sliderSteps = makeSlider(ages, agesMap);
    traces = makeTraces(locs, ordAges, locsColors, cparameters);
    
    var myPlot = document.getElementById('graph'),
        data = traces["mrk"].concat(traces["time"])
        layout = {margin: {l: 5, r: 5, t: 5, b: 5},
                  paper_bgcolor: "#e7e7e7",
                  width: fdims[0],
                  height: fdims[1],
                  hovermode:'closest',
                  grid: {rows: 2, columns: 1, pattern: 'independent'},
                  xaxis2: {
                      anchor: 'y2', domain: [0, 1], range: [tticks[0], tticks[tticks.length-1]],
                      tickmode: "array", tickvals: tticks, ticktext: ttickTxts,
                      showgrid: true, showline: false, zeroline: false, showticklabels: true
                  },
                  yaxis2: {
                      anchor: 'x2', domain: [0, 0.44], // range: tlims,
                      showgrid: false, showline: false, zeroline: false, showticklabels: false
                  },
                  geo: {
                      scope: 'world',
                      resolution: 50,
                      projection: {type: 'miller'},
                      domain: {x: [0, 1], y: [0.46, 1]},
                      lonaxis: {range: xlims},
                      lataxis: {range: ylims},
                      showcoastlines: true,
                      showocean: true,
                      oceancolor: "#F9FCFF",
                      showcountries: true
                  },
                  // xaxis: {
                  //     showgrid: true,
                  //     gridwidth: 0.5,
                  //     range: xlims,
                  //     dtick: 5
                  // },
                  // yaxis: {
                  //     showgrid: true,
                  //     gridwidth: 0.5,
                  //     range: ylims,
                  //     dtick: 5
                  // },
                  sliders: [{
                      pad: {l: 0, t: 40},
                      currentvalue: {
                          visible: true,
                          prefix: 'Time (Ma): ',
                          xanchor: 'left',
                          font: {size: 20, color: '#6b6b6b'}
                      },
                      steps: sliderSteps
                      }]
                  };
    config = {displaylogo: false,
              responsive: true,
              toImageButtonOptions: {
                  format: 'svg', // one of png, svg, jpeg, webp
                  filename: 'mapFossils',
                  width: fdims[0],
                  height: fdims[1],
                  scale: 1 // Multiply title/legend/axis/canvas sizes by this factor
              },
              modeBarButtonsToRemove:  ['zoom2d', 'pan2d', 'select2d', 'lasso2d', 'zoomIn2d', 'zoomOut2d',
                                        'autoScale2d', 'resetScale2d',
                                        'zoomInGeo', 'zoomOutGeo', 'resetGeo', 'hoverClosestGeo',
                                        'hoverClosestGl2d', 'hoverClosestPie', 'toggleHover', 'resetViews',
                                        'sendDataToCloud', 'toggleSpikelines', 'resetViewMapbox',
                                        'hoverClosestCartesian', 'hoverCompareCartesian']
             }
    
    Plotly.newPlot('graph', data, layout, config);

    // function highlightOn(points, nbRows){
    //     for(var i=0; i < points.length; i++){
    //         bidsOn = getAssociatedIds(points[i].curveNumber, nbRows);
    //         if ( bidsOn.length == 2 ){
    //             Plotly.restyle('graph', {'marker.size': [2*baseMrkSz]}, [bidsOn[0]]);
    //             Plotly.restyle('graph', {'line.width': [2*baseLineW]}, [bidsOn[1]]);
    //         }
    //     }
    // }
    // function highlightOff(points, nbRows){
    //     for(var i=0; i < points.length; i++){
    //         bidsOff = getAssociatedIds(points[i].curveNumber, nbRows);
    //         if ( bidsOff.length == 2 ){
    //             Plotly.restyle('graph', {'marker.size': [baseMrkSz]}, [bidsOff[0]]);
    //             Plotly.restyle('graph', {'line.width': [baseLineW]}, [bidsOff[1]]);
    //         }            
    //     }
    // }
    // myPlot.on('plotly_hover', function(data){ highlightOn(data.points, locs.lidnum.length); });
    // myPlot.on('plotly_unhover', function(data){ highlightOff(data.points, locs.lidnum.length); });
}

function loadGraph(cparameters){
    // console.log(cparameters);
    
    d3.csv(cparameters["locfile"], function (d) { return rtrnLoc(d, cparameters["v_count"], cparameters["v_focus"]) })
    .then(
        function(dataLocs){ makePlot(dataLocs, cparameters); },
        function(err) { console.log(err); } // Error: "It broke when reading supps"
    );

}




regionLbls = ["Europe-wide", "North-America"];
regionIdsMap = {"Europe-wide": "EU", "North-America": "NA"};
faunaLbls = ["Large", "Small", "Carnivorous"];
faunaIdsMap = {"Large": "L", "Small": "S", "Carnivorous": "C"};
// prepare drop-down selectors
var regionSelector = document.querySelector('#choice-region');
assignOptions(regionLbls, regionSelector);
var faunaSelector = document.querySelector('#choice-fauna');
assignOptions(faunaLbls, faunaSelector);
var focusSelector = document.querySelector('#choice-focus');

function updateDataSrcRegion(){
    updateDataSrc(1, 0, 0);
}
function updateDataSrcFauna(){
    updateDataSrc(0, 1, 0);
}
function updateDataSrcFocus(){
    updateDataSrc(0, 0, 1);
}
function updateDataSrc(upRegion, upFauna, upFocus){

    var queryString = window.location.search;
    var urlParams = new URLSearchParams(queryString);
    if ( urlParams.has('region') && regionLbls.includes(urlParams.get('region')) ) {
        regionSelector.value = urlParams.get('region');
    }
    if ( urlParams.has('fauna') && faunaLbls.includes(urlParams.get('fauna')) ) {
        faunaSelector.value = urlParams.get('fauna');
    }
    
    var regionId = regionIdsMap[regionSelector.value];
    var faunaId = faunaIdsMap[faunaSelector.value];
    
    var focusLbls = ["MEAN_HYPSODONTY"].concat(foci[regionId+"-"+faunaId]);
    if ( upRegion + upFauna > 0){
        clearOptions(focusSelector);
        assignOptions(focusLbls, focusSelector);
    }
    if ( urlParams.has('vfocus') && focusLbls.includes(urlParams.get('vfocus')) ) {
        focusSelector.value = urlParams.get('vfocus');
    }
    if ( urlParams.has('vfocus') && focusLbls.includes(tmpFocus) ) {
        focusSelector.value = tmpFocus;
    }

    var dfocus = "MEAN_HYPSODONTY";   
    var vfocus = focusSelector.value;
    var vmin = 1;
    var vmax = 3;    
    if ( vfocus != dfocus ) {
        vmax = -1;
    }
    // console.log(vfocus, vmax);
          
    var vcount = "nb_genera";
    var cmap = "summer";
    if ( urlParams.has('vcount') ) {
        vcount = urlParams.get('vcount');
    }
    if ( urlParams.has('vmax') ) {
        vmax = parseInt(urlParams.get('vmax'));
    }
    if ( urlParams.has('vmin') ) {
        vmin = parseInt(urlParams.get('vmin'));
    }
    if ( urlParams.has('cmap') ) {
        cmap = urlParams.get('cmap');
    }
    
    var loc_parameters = {"fdims": [1200, 800],
                          "baseMrkSz": 3,
                          "noiseLoc": 0,
                          "padT": 0.001                          
                         }

    loc_parameters["cmap"] = cmap;
    
    loc_parameters["region_id"] = regionId;
    loc_parameters["fauna_id"] = faunaId;
    loc_parameters["v_focus"] = vfocus;

    loc_parameters["v_min"] = vmin;
    loc_parameters["v_max"] = vmax;
    
    loc_parameters["v_count"] = vcount;
    loc_parameters["locfile"] = "data/outline_"+regionId+"-"+faunaId+".csv";
    if ( regionId == "EU" ) {
        loc_parameters["xlims"] = [ -15, 76];
        loc_parameters["ylims"] = [ 30, 60];
        loc_parameters["slice_names"] = ["MN2", "MN3", "MN4","MN5", "MN6","MN7-8", "MN9", "MN10","MN11","MN12","MN13", "MN14","MN15", "MN16"];
        loc_parameters["ttickTxts"] = ["","MN2","21","MN3", "18", "MN4","MN5", "15", "MN6","MN7-8","12", "MN9", "MN10","MN11","MN12","MN13", "6","MN14","MN15", "MN16", "3", ""];
        loc_parameters["agesTicks"] = [21.7, 21, 19.5, 18, 17.2, 16.4, 15, 14.2,12.85, 12, 11.2, 9.9,8.9,7.6,7.1, 6, 5.3, 5, 3.55, 3];
        loc_parameters["baseLineW"] = .8;
    } else {
        loc_parameters["xlims"] = [ -130, -65];
        loc_parameters["ylims"] = [ 21, 86];
        loc_parameters["slice_names"] = ["Orellan","Whitneyan", "Arikareean-1","Arikareean-2", "Arikareean-3","Arikareean-4","Hemingfordian-1","Hemingfordian-2","Barstovian-1","Barstovian-2", "Clarendonian-1","Clarendonian-2","Clarendonian-3","Hemphillian-1","Hemphillian-2","Hemphillian-3","Hemphillian-4","Blancan-Early"];
        loc_parameters["ttickTxts"] = ["","Orellan","Whitneyan", "31","Arikareean-1","Arikareean-2", "25","Arikareean-3","Arikareean-4","Hemingfordian-1","Hemingfordian-2","Barstovian-1","Barstovian-2", "13", "Clarendonian-1","Clarendonian-2","Clarendonian-3","Hemphillian-1","Hemphillian-2","Hemphillian-3","Hemphillian-4","Blancan-Early","3", ""];
        loc_parameters["agesTicks"] = [33.9,32.1, 31,29.75,28, 25,22.2,19.09,18.5,17.34,15.97,14.78,13, 12.57,12.11,10.09,9.07,7.59,6.88,5.91,4.91, 3];
        loc_parameters["baseLineW"] = 1.5;
    }
    loadGraph(loc_parameters);
}
regionSelector.addEventListener('change', updateDataSrcRegion, false);
faunaSelector.addEventListener('change', updateDataSrcFauna, false);
focusSelector.addEventListener('change', updateDataSrcFocus, false);

updateDataSrc(1,1,1);
