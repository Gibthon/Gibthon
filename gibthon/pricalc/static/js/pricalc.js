// Validation for input

function validate_gene(gene){
	var error = gene.nextElementSibling.nextElementSibling;
	var genec = error.nextElementSibling;
	
	var dnare = new RegExp("^[atgc]*$", 'i');

	
	if(!dnare.test(gene.value)){
		writeerror(error, "Sequence must only contain atgc (case insensitive)", 1);
		validate(gene.form);
		return false;
	}
	if(gene.value.length < 20){
		writeerror(error, "Sequence too short, must be at least 20 characters", 1);
		validate(gene.form);
		return false;
	}
	if(gene.value.length > 100){
		writeerror(error, "Sequence too long, must be at most 100 characters", 1);
		validate(gene.form);
		return false
	}
	writeerror(error, "OK!", -1);
	genec.value = gene.value;
	validate(gene.form);
	return true;
}

function writeerror(error, errortext, status){
	error.innerHTML = errortext;
	if(status == 1){
		error.setAttribute("class", "error");
		$("#tabs").tabs("option", "disabled", [1]);
	}else if (status == -1){
		error.setAttribute("class", "errorok");
		validate(error.parentNode.form);
	}else if (status == 0){
		error.setAttribute("class", "errorwait");
		$("#tabs").tabs("option", "disabled", [1]);
	}
}



function validate(form){
	with(form){	
		if(lgene.nextElementSibling.nextElementSibling.getAttribute("class") == "errorok" && rgene.nextElementSibling.nextElementSibling.getAttribute("class") == "errorok"){
			$("#tabs").tabs("option", "disabled", []);
			return true;
		}
	}
	return false;
}

// select everything in the box
function selectall(event){
	var text_val=event.currentTarget;
	text_val.focus();
	text_val.select();
}

var minlength = 20;

// segment of oligo
var Segment = function (side, gene, start, end){
	// left or right side
	this.side = side;
	// gene sequence of whole thing
	this.gene = gene;
	// start of oligo section
	this.start = start;
	// end of oligo section
	this.end = end;
	// melting temperature
	this.Tm = null;
}

// for printing the interface to gene
Segment.prototype.print_gene = function(){
	// this is where we're putting the gene
	var out = document.getElementById(this.side + 'gene');
	// break it down now y'all
	var g = this.gene.split("");
	// lets add everything into one variable, and then put it in in one go
	var outs = '<ul class="gene" id="' + this.side + '">';
	for (var i = 0; i < g.length; i++){
		// counting up or down?
		var n = (this.side == 'l') ? g.length - i : i+1;
		if(n == this.start){
			// start base
			outs = outs + liHelper("start droppable start" + this.side, g[i], n);
		}else if (n == this.end){
			// end base
			outs = outs + liHelper("end droppable end" + this.side, g[i], n);
		}else if (n > this.end - minlength + 1 && n < this.start + minlength - 1){
			// forbidden base
			outs = outs + liHelper("invalid droppable", g[i], n);
		}else if (n < this.end && n > this.start){
			// base in gene
			outs = outs + liHelper("middle droppable", g[i], n);
		}else{
			// boring old base
			outs = outs + liHelper("droppable", g[i], n);
		}
	}
	// close list
	outs = outs + '<\/ul>';
	// gogogo!
	out.innerHTML=outs;	
}

Segment.prototype.recut = function(index, from){
	if(from=="start"){
		if(index > this.end - minlength + 1){
			move_error(this.side, from, index);
			return false;
		}else{
			this.start = index;
		}
	}else{
		if(index < this.start + minlength - 1){
			move_error(this.side, from, index);
		}else{
			this.end = index;
		}
	}
}
		
		
// for printing each <li> bit of printgene - saves repetition
function liHelper(classes, g, n){
	return '<li class="' + classes + '"><a href="#" onclick="return false;" class="' + n + '">' + g + '<\/a><\/li>';
}

// error for making gene too short!
function move_error(side, end, index){
	$('#prompt').prompt({
		title:'Warning',
		message:"Can't move " + end + " of " + (side == 'l' ? "left" : "right") + " gene to position " + index + ": overlap too short"
	});
}

// oligo object
var Oligo = function (leftgene, rightgene){
	this.L = new Segment('l', leftgene, 1, 20);
	this.R = new Segment('r', rightgene, 1, 20);
	this.Tm = null;
	this.dG = null;
	this.lseq = null;
	this.rseq = null;
	this.id = null;
	this.alwaysdg = document.getElementById('alwaysdg');
	this.mg = document.getElementById('mgsalt');
	this.na = document.getElementById('nasalt');
	this.get_info();
}

Oligo.prototype.get_info = function(){
	$.ajax({
		type:'post',
		url:'go',
		data:{
			gene1:this.L.gene,
			gene2:this.R.gene,
			el:this.L.end,
			sl:this.L.start,
			er:this.R.end,
			sr:this.R.start
		},
		context:this,
		dataType:'json',
		success: function(data) {
			this.L.Tm = data.TmT;
			this.R.Tm = data.TmB;
			this.Tm = data.TmF;
			this.lseq = data.SeqT;
			this.rseq = data.SeqB;
			$('#time').html((data.time*1000).toString().substring(0,4) + "ms");
			this.print_info();
		}
	});
}

// http_request maker
// creates the http function for the oligo, and sets up the return handler 
//  - ie, dumps data into oligo and calls print function
/*function getH(parent) {
	var h = new XMLHttpRequest();
	h.parent = parent;
	h.onreadystatechange = function() {
	  	if (this.readyState == 4){
			var JSON_ret = JSON.parse(this.responseText);
			if (JSON_ret['error'] < 0){
				this.parent.L.Tm = JSON_ret['Tml'];
				this.parent.R.Tm = JSON_ret['Tmr'];
				this.parent.Tm = JSON_ret['TmO'];
				this.parent.dG = JSON_ret['dG'];
				this.parent.lseq = JSON_ret['lseq'];
				this.parent.rseq = JSON_ret['rseq'];
				this.parent.id = JSON_ret['id'];
				this.parent.print_info();
			}
		}
	}
	return h;
}*/

Oligo.prototype.print_info = function(){
	this.L.print_gene();
	this.R.print_gene();
	document.getElementById('tmleft').innerHTML=this.L.Tm.toString().substring(0,4) + "&deg;C";
	document.getElementById('tmright').innerHTML=this.R.Tm.toString().substring(0,4) + "&deg;C";
	document.getElementById('tmall').innerHTML=this.Tm.toString().substring(0,4) + "&deg;C";	
	document.getElementById('lseq').innerHTML=this.lseq;
	document.getElementById('rseq').innerHTML=this.rseq;
	var dgb = document.getElementById('dg');
	if(!this.dG){
		if(dgb.innerHTML.indexOf("Calculate") < 0 ){
			if(dgb.innerHTML.indexOf("Update") < 0 && dgb.innerHTML.indexOf('.gif') < 0){
				var old = dgb.innerHTML;
				dgb.innerHTML = old + ' <a href="#" onClick="MyOligo.get_dg(); return false">Update</a>';
			}
		}
	}else{
		document.getElementById('pic').innerHTML="<img src='out/" + this.id + ".jpg' />"
		dgb.innerHTML = this.dG;
	}
	
	var that = this;
	$(".end").draggable({
        helper: function(event) {
            return $('<div class="endhelp"></div>');
        },
        grid:[20,30],
        containment: 'parent'
    });
	$(".start").draggable({
        helper: function(event) {
            return $('<div class="starthelp"></div>');
        }, 
        grid:[20,30], 
        containment:'parent'
    });
	$(".droppable").droppable({
		drop: function(event, ui) {
			var gene = event.target.parentNode.id;
			var index = event.target.firstChild.className;
			var from = ui.draggable.context.className.split(" ")[0];
			that.recut(gene, parseInt(index), from);
		}
	});
}

Oligo.prototype.recut = function(gene, index, from){
	if(gene=='l'){
		this.L.recut(index,from);
	}else{
		this.R.recut(index,from);
	}
	if(this.alwaysdg.checked){
		this.get_dg();
	}else{
		this.get_info();
	}
}

Oligo.prototype.get_dg = function(){
//	this.http.open("GET", "cgi-bin/gibthon.cgi?action=godg&gene1=" + this.L.gene + "&gene2=" + this.R.gene + "&el="+this.L.end+"&sl="+this.L.start+"&er="+this.R.end+"&sr="+this.R.start+"&mg="+this.mg.value+"&na="+this.na.value,true);
//	this.http.send(null);
	$('#dg').html('<img src="'+STATIC_URL+'/images/spinner.gif"/>');
}


function saltChange(){
	var dgb = document.getElementById('dg');
	if(dgb.innerHTML.indexOf("Calculate") < 0 ){
		if(dgb.innerHTML.indexOf("Update") < 0 && dgb.innerHTML.indexOf('.gif') < 0){
			var old = dgb.innerHTML;
			dgb.innerHTML = old + ' <a href="#" onClick="MyOligo.get_dg(); return false">Update</a>';
		}
	}
}

function printtime(t){
	$('#time').html((t*1000).toString().substring(0,6) + "ms");
}

function getcsv(){
	var dgb = document.getElementById('dg');
	if(dgb.innerHTML.indexOf("Calculate") > 0 || dgb.innerHTML.indexOf("Update") > 0){
		$('#prompt').prompt({
			title:'Note',
			message:"Please refresh dg first"
		});
		return false;
	}
	window.open("csv.php?seq=" + JSON_ret['full'] + "&tl=" + JSON_ret['Tml'].toString().substring(0,4) + "&tr=" + JSON_ret['Tmr'].toString().substring(0,4) + "&t=" + JSON_ret['TmO'].toString().substring(0,4) + "&dg=" + JSON_ret['dG']);
	return false;
}