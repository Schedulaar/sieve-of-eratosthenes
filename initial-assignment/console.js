var groupBy = function(xs, key) {
  return xs.reduce(function(rv, x) {
    (rv[x[key]] = rv[x[key]] || []).push(x);
    return rv;
  }, {});
};

const pData = groupBy(data, 'p')

function objectMap(object, mapFn) {
  return Object.keys(object).reduce(function(result, key) {
    result[key] = mapFn(object[key])
    return result
  }, {})
}

const procNData = objectMap(pData, val => groupBy(val, 'n'))

const meansAndDeviations = objectMap(procNData, pdata => objectMap(pdata, iterations => {
  const mean = iterations.reduce((acc,val) => acc + val.t,0)/iterations.length
  const sum = iterations.reduce((acc,val) => acc + (val.t - mean)*(val.t - mean),0)
  const relDeviation = Math.sqrt(sum / (iterations.length - 1)) / mean
  return { mean, relDeviation }
}))

const latexData = objectMap(meanAndDeviationObj2, pdata => {
  const means = Object.keys(pdata).reduce((acc, n) => acc += `(${n}, ${pdata[n].mean})
`, '')
  const relDeviations = Object.keys(pdata).reduce((acc, n) => acc += `(${n}, ${pdata[n].relDeviation})
`, '')
  return { means, relDeviations }
})

const means = Object.values(latexData).reduce((acc, val) => acc + `

\\addplot coordinates {
${val.means}
};`, '')

const relDeviations = Object.values(latexData).reduce((acc, val) => acc + `

\\addplot coordinates {
${val.relDeviations}
};`, '')