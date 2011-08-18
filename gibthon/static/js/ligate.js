$(document).ready(function() {
	react = new Reaction();
	$(":input").change(function(){
		react.go();
	});
	$('#button-addrow').button({	icons:{primary:'ui-icon-plusthick'},
						label: "Add Row"
					});
	$('#button-return').button({ icons:{secondary:'ui-icon-home'}, label:'Home', text:false});
});



var Reaction = function(){
	var forms = document.getElementById('input').children[0].children;
	this.components = [];
	for (x in forms){
		if (typeof forms[x]==='object'){
			if(x!='0'){
				this.components.push(new Component(forms[x]));
			}
		}
	}
	var setform = document.getElementById('base');
	this.tv = setform.total_volume;
    this.bv = setform.buffer_volume;
    this.lv = setform.ligase_volume;
    this.maxvol = '';
    this.go();
}

Reaction.prototype.validate = function(){
	var numre = /^\d+(?:\.\d*)?$/;
	var ret = true;
	for (var i = 0; i < this.components.length; i++){
		if (!this.components[i].validate()){
			ret = false;
		}
	}
	ret = ret*valcheck(numre, this.tv);
	ret = ret*valcheck(numre, this.bv);
	ret = ret*valcheck(numre, this.lv);
	if(!ret){
		return ret;
	}
	this.maxvol = parseFloat(this.tv.value) - parseFloat(this.bv.value) - parseFloat(this.lv.value);
	if (parseFloat(this.maxvol) <= 0){
		alert("Total volume too small!");
		return false;
	}
	return ret;
}

Array.prototype.max=function()
{
	if (this.length == 0)
		return {'index':-1};
	var maxIndex = 0;
	for (var i = 1; i<this.length; i++)
		if (this[i]>this[maxIndex]) 
			maxIndex = i;
	return {'index':maxIndex, 'value':this[maxIndex]};
}


Reaction.prototype.go = function(){
	if (!this.validate()){
		return false;
	}
	var comps = [];
	var tmol = 0;
	var doms = [];
	
	for (var i = 0; i < this.components.length; i++){
		tmol += this.components[i].mol;
		doms.push(this.components[i].mol/parseFloat(this.components[i].maxvol.value));
		comps.push(i);
	}
	var dom = doms.max().index;
	this.components[dom].vol = Math.min(parseFloat(this.components[dom].maxvol.value), (this.maxvol * this.components[dom].mol)/tmol);
	comps.splice(dom,1);
	for (var i = 0; i < comps.length; i++){
		this.components[comps[i]].vol = (this.components[dom].vol / this.components[dom].mol) * this.components[comps[i]].mol;
	}
	out = document.getElementById('output');
	outs = "<table style=\"width:125px;\"><tr><th style=\"width:100px;\"></th><th style=\"width:25px;\">Vol (µl)</th></tr><tr><td>Buffer</td><td>" + this.bv.value + "</td></tr><tr><td>Enzyme</td><td>" + this.lv.value + " </td></tr>";
	var water = this.maxvol;
	for (var i = 0; i < this.components.length; i++){
		outs += "<tr><td>"+this.components[i].name.value+"</td><td>"+round2dp(this.components[i].vol)+"</td></tr>";
		water -= this.components[i].vol;
	}
	outs += "<tr><td>Water</td><td>"+round2dp(water)+"</td></tr><tr><th>Total</th><td>"+this.tv.value+"</td></tr></table>";
    out.innerHTML = outs;
}

function round2dp(val){
	return Math.round(val*10)/10;
}

var Component = function(form){
	with(form){
		this.name = children[0].children[0];
		this.ratio = children[1].children[0];
		this.conc = children[2].children[0];
		this.length = children[3].children[0];
		this.maxvol = children[4].children[0];
	}
	this.uperul = '';
	this.mol = '';
	this.vol = '';
}

Component.prototype.validate = function(){
	var numre = /^\d+(?:\.\d*)?$/;
	var intre = /^\d+$/;
	var ret = true;
	ret=ret*valcheck(numre, this.ratio);
	ret=ret*valcheck(numre, this.conc);
	ret=ret*valcheck(intre, this.length);
	ret=ret*valcheck(numre, this.maxvol);
	if(!ret){
		return ret;
	}
	this.uperul = parseFloat(this.conc.value)/parseFloat(this.length.value);
	this.mol = parseFloat(this.ratio.value)/this.uperul;
	return ret;
}

function addTableRow(jQtable){
    jQtable.each(function(){
        var $table = $(this);
        var tds='<tr><td><input type="text" name="name" value="Insert" />\
        </td><td><input type="text" name="ratio" value="3" class="numeric" />\
        </td><td><input type="text" name="conc" class="numeric" value="20"/>\
        </td><td><input type="text" name="length" class="numeric" value="1000"/>\
        </td><td><input type="text" name="maxvol" class="numeric" value="10" />\
        </td><td><a href="#" onClick="rmTableRow($(this)); return false;"><span class="ui-icon ui-icon-trash"></span></a></td></tr>'
        if($('tbody', this).length > 0){
            $('tbody', this).append(tds);
        }else {
            $(this).append(tds);
        }
    });
    react = new Reaction();
    handleAutoComplete();
    $(":input").change(function(){
		react.go();
	});
}

function rmTableRow(row){
	$(row[0].parentNode.parentNode).remove();
	react = new Reaction();
}

function valcheck(re, id){
	if(re.exec(id.value) === null){
		id.parentNode.className="errorcell";
		return false;
	}else{
		id.parentNode.className="";
		return true;
	}
}
