let amplifier=40;
let patche0 =[[0,0.262866,0.4253250],
          [0, -0.262866, 0.4253250],
          [0.425325, 0, 0.2628660],
          [0, -0.262866, -0.4253250],
          [0, 0.262866, -0.4253250],
          [-0.425325, 0, -0.2628660]]
let colors0=[21,22,23,24,25,26];
let colors=[[255,102,102],[51,119,255],[204,204,0],[43,255,0],[204,102,255],[51,255,221]]
let m=[]
let i=0;
patche0.forEach(vec => {
  vec = vec.map(function(x){return x*amplifier});
  //console.log(vec);
  let p = new THREE.Vector3(vec[0],vec[1],vec[2]);
    m[i] = new THREE.Mesh(
    new THREE.SphereGeometry(5,64,64),
    new THREE.MeshPhongMaterial(
      {color: new THREE.Color(colors[i][0]/256.0,colors[i][1]/256.0,colors[i][2]/256.0),
        opacity:.5
      }
    )
  )
  //m.material.transparent=true;
  m[i].position.copy(p);
  scene.add(m[i]);
  i+=1;
});
