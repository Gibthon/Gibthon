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

var minlength = 10;

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
	var out = $('#'+this.side + 'gene');
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
	out.html(outs);	
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
	this.lseq = null;
	this.rseq = null;
	this.id = null;
	this.mg = $('#mgsalt')[0];
	this.na = $('#nasalt')[0];
	this.dgbutton = $('button#self-prime-check-button').button();
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
			sr:this.R.start,
			mg:this.mg.value,
			na:this.na.value
		},
		context:this,
		dataType:'json',
		success: function(data) {
			this.L.Tm = data.TmT;
			this.R.Tm = data.TmB;
			this.Tm = data.TmF;
			this.lseq = data.SeqT;
			this.rseq = data.SeqB;
			this.print_info();
		}
	});
}

Oligo.prototype.get_ss = function() {
	$.ajax({
		type:'post',
		url:'selfprime',
		data:{
			gene1:this.L.gene,
			gene2:this.R.gene,
			el:this.L.end,
			sl:this.L.start,
			er:this.R.end,
			sr:this.R.start,
			mg:this.mg.value,
			na:this.na.value
		},
		context:this,
		dataType:'json',
		success: function(data) {
			$('#boxplot').html('<img src="'+data.image+'" />');
			this.dgbutton.button('enable');
		}
	});
}
	

Oligo.prototype.print_info = function(){
	this.L.print_gene();
	this.R.print_gene();
	$('#tmleft').html(this.L.Tm.toString().substring(0,4) + "&deg;C");
	$('#tmright').html(this.R.Tm.toString().substring(0,4) + "&deg;C");
	$('#tmall').html(this.Tm.toString().substring(0,4) + "&deg;C");
	$('#lseq').html(this.lseq);
	$('#rseq').html(this.rseq);
	$('#length').html(this.lseq.length);
	
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
	this.get_info();
}


function saltChange(){
	var dgb = $('#dg');
	if(dgb.html().indexOf("Calculate") < 0 ){
		if(dgb.html().indexOf("Update") < 0 && dgb.html().indexOf('.gif') < 0){
			dgb.html(dgb.html() + ' <a href="#" onClick="MyOligo.get_dg(); return false">Update</a>');
		}
	}
}

$(document).ready(function () {
	$('#tabs').tabs({
		disabled:[1],
		show: function(event, ui) {
			if(ui.index == 1){
				var form = $('#form')[0];
				MyOligo = new Oligo(form.lgenec.value.toLowerCase(), form.rgenec.value.toUpperCase());
			}
		}
	});

	 $('button#self-prime-check-button').click(function() {
		$(this).button('disable');
		MyOligo.get_ss();
	});

	$('input[name="lgene-search"]').autocomplete({
		source:'list',
		minLength:1,
		focus: function( event, ui ){
			return false;
		},
		select:function( event, ui ){
			$('textarea[name="lgene"]')[0].value = ui.item.last.toLowerCase();
			validate_gene($('textarea[name="lgene"]')[0]);
		}
	});
	$('input[name="rgene-search"]').autocomplete({
		source:'list',
		minLength:1,
		focus: function( event, ui ){
			return false;
		},
		select:function( event, ui ){
			$('textarea[name="rgene"]')[0].value = ui.item.first.toLowerCase();
			validate_gene($('textarea[name="rgene"]')[0]);
		}
	});
	
	$('#left-up').click(function () {
		if (MyOligo.L.end < $('#lgene ul').children().length) {
			MyOligo.L.end += 1;
			MyOligo.get_info();
		}
	});
	$('#left-down').click(function () {
		if (MyOligo.L.end - MyOligo.L.start > $('#lgene ul li[class^="invalid"]').length + 2 ){
			MyOligo.L.end -= 1;
			MyOligo.get_info();
		}
	});
	$('#right-up').click(function () {
		if (MyOligo.R.end < $('#rgene ul').children().length) {
			MyOligo.R.end += 1;
			MyOligo.get_info();
		}
	});
	$('#right-down').click(function () {
		if (MyOligo.R.end - MyOligo.R.start > $('#rgene ul li[class^="invalid"]').length + 2 ){
			MyOligo.R.end -= 1;
			MyOligo.get_info();
		}
	});
	$('.all-up').click(function () {
		if (MyOligo.L.end < $('#lgene ul').children().length && MyOligo.R.end < $('#rgene ul').children().length) {
			MyOligo.R.end += 1;
			MyOligo.L.end += 1;
			MyOligo.get_info();
		}
	});
	$('.all-down').click(function() {
		if (MyOligo.R.end - MyOligo.R.start > $('#rgene ul li[class^="invalid"]').length + 2  && MyOligo.L.end - MyOligo.L.start > $('#lgene ul li[class^="invalid"]').length + 2) {
			MyOligo.R.end -= 1;
			MyOligo.L.end -= 1;
			MyOligo.get_info();
		}
	});
});

