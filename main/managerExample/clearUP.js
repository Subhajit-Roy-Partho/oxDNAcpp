// select polyTS then their 3' ends. After that drag drop this.

function deleteIds(ids) { // delete based on ids
  let ele = api.getElements(ids);
  edit.deleteElements(ele);
}

let juncturePoints =[];
let i=0;
selectedBases.forEach(base => {// Do the operations on 3' selected bases
  let id = base.id;
  id -=24;
  id=api.getElements([id])[0].n5.id;
  juncturePoints.push(id);

  deleteIds([...Array(25).keys()].map(x => x+base.id-24));

});

deleteIds([...Array(21).keys()].map(x => x+15930))//extra strand no use

api.selectElementsIDs(juncturePoints);
console.log(juncturePoints);
render();

// 12882, 13397,13366
// 13130, 13925, 12582