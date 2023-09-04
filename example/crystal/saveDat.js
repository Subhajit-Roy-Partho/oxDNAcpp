function prepareDat() {
  let dat="";
  dat=[
    't = 0',
    `b = ${Math.ceil(box.x)} ${Math.ceil(box.x)} ${Math.ceil(box.x)}`,
    'E = 0 0 0\n'
  ].join('\n');
  elements.forEach(e => {
    dat+=e.getDatFileOutput();
  });
  return dat;
}

function download(filename, text) {
    var element = document.createElement('a');
    element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
    element.setAttribute('download', filename);

    element.style.display = 'none';
    document.body.appendChild(element);

    element.click();

    document.body.removeChild(element);
}

function repulsion_plane(position, dir){
  return `\n{\ntype = repulsion_plane\nparticle = -1\nstiff = 1\ndir = ${dir[0]}, ${dir[0]}, ${dir[0]}\nposition = ${position}\n}`
}

function trap(particle, position){
  return `\n{\ntype = trap\nparticle = ${particle}\npos0 = ${position.x}, ${position.y}, ${position.z}\nstiff = 1\nrate =0\ndir = 1,0,0\n}`
}

function forceDat(){
  let force="";
  force+=repulsion_plane(0,[0,0,1]);
  force+=repulsion_plane(Math.ceil(box.x),[0,0,-1]);
  selectedBases.forEach(e=>{
    force+=trap(e.id,e.getPos());
  })
  return force;
}


download('test.txt',prepareDat());
download("forces.txt", forceDat());