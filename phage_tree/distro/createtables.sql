drop database phage; create database phage; use phage;
create table phage (
	count INT NOT NULL PRIMARY KEY AUTO_INCREMENT,
	organism VARCHAR(255) NOT NULL,
	genbankname VARCHAR(255) NOT NULL,
	accession VARCHAR(255) NOT NULL,
	locus VARCHAR(255) NOT NULL,
	beginning INT,
	end INT,
	sequence LONGTEXT,
	family VARCHAR(255),
	phylogeny LONGTEXT
	);

create table protein (
	count INT NOT NULL PRIMARY KEY AUTO_INCREMENT,
	organism INT NOT NULL,
	start INT,
	stop INT,
	complement INT,
	gene VARCHAR(255),
	function  VARCHAR(255),
	note LONGTEXT,
	product VARCHAR(255),
	proteinid VARCHAR(255),
	dbxref VARCHAR(255),
	translation LONGTEXT
	);

create table swiss (
	count int PRIMARY KEY NOT NULL auto_increment,
	ac longtext NULL,
	cc longtext NULL,
	de longtext NULL,
	dr longtext NULL,
	dt longtext NULL,
	ft longtext NULL,
 	gn longtext NULL,
 	id longtext NULL,
 	kw longtext NULL,
 	oc longtext NULL,
 	og longtext NULL,
 	os longtext NULL,
 	ra longtext NULL,
 	rc longtext NULL,
 	rl longtext NULL,
 	rn longtext NULL,
 	rp longtext NULL,
 	rt longtext NULL,
 	rx longtext NULL,
 	sq longtext NULL,
 	seq longtext NULL
);
grant all on phage to rob;
grant all on protein to rob;
grant all on swiss to rob;
