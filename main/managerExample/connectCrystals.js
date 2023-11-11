handle1 = new THREE.Vector3(0,0,0);
handle1Elements = api.getElements([12882, 13397,13366]);
// handle2Elements = api.getElements([13130, 13925, 12582]);

for(let i=0;i<handle1Elements.length;i++){
    handle1.add(handle1Elements[i].getPos());
}
handle1.divideScalar(handle1Elements.length);

centeroid = systems[0].getCom();

unitVector = centeroid.sub(handle1);

systems[0].select;

cutWrapper();

pasteWrapper(true);

selectionToCluster();

// centeroid = systems[0].getCom();

// centroid.add(THREE.Vector3(80,80,80));
// translateElements(selectedBases,centeroid);