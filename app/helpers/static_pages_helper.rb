module StaticPagesHelper

	require 'helper.rb'
	require 'proteinFamily.rb'

	Anamorphs = {"Cll" => "Candida lusitaniae", "Cyj" => "Candida utilis",
		"Deh" => "Candida famata var. flareri", "Deb" => "Brettanomyces bruxellensis",
		"Ke" => "Candida holmii", "Kl" => "Candida sphaerica",
		"Km" => "Candida kefyr", "Mrg" => "Candida guilliermondii",
		"Piu" => "Candida acidothermophilum", "Phm" => "Candida valida",
		"Sc_an" => "Candida robusta", "Sc_aq" => "Candida robusta",
		"Sc_al" => "Candida robusta", "Sc_b" => "Candida robusta", 
		"Sc_a" => "Candida robusta", "Sc_c" => "Candida robusta",
		"Sc_am" => "Candida robusta", "Tod" => "Candida colliculosa",
		"Yl" => "Candida lipolytica"
	}

	Altnames = {"Bai" => "Pichia inositovora, Yamadazyma inositovora", 
		"Ca_a" => "Candida claussenii (obsolete), Candida langeroni (obsolete), Candida stellatoidea (obsolete), Monilia albicans (obsolete)",
		"Ca_b" => "Candida claussenii (obsolete), Candida langeroni (obsolete), Candida stellatoidea (obsolete), Monilia albicans (obsolete)",
		"Cgl" => "Cryptococcus glabratus (obsolete), Torulopsis glabrata (obsolete), Torulopsis stercoralis (obsolete)",
		"Cao" => "Candida parapsilosis group II",
		"Ct_a" => "Candida paratropicalis (obsolete), Candida vulgaris (obsolete), Monilia tropicalis (obsolete)",
		"Cll" => "Clavispora imtechensis (misnomer)",
		"Cyj" => "Pichia jadinii, Candida guilliermondii var. nitratophila, Hansenula jadinii, Torula utilis, Torulopsis utilis, Pichia jadinii (Candida utilis), Lindnera jadinii (misnomer)",
		"Deb" => "Brettanomyces custersii",
		"Erg" => "Ashbya gossypii, Eremothecium gossypii (Ashby et Nowell) Kurtzman, Nematospora gossypii",
		"Hyb" => "Endomycopsis burtonii",
		"Kaa" => "Kluyveromyces africanus",
		"Ke" => "Saccharomyces exiguus",
		"Ka" => "Kluyveromyces aestuarii MY 111",
		"Kl" => "Kluyveromyces marxianus var. lactis, Kluyveromyces marxianus var. drosophilarum, Kluyveromyces marxianus lactis, Kluyveromyces drosophilarum, Kluyveromyces lactis var. drosophilarum, Kluyveromyces lactis var. lactis",
		"Km" => "Kluyveromyces marxianus var. marxianus, Kluyveromyces fragilis, Kluyveromyces cicerisporus, Candida pseudotropicalis",
		"Kop_a" => "Pichia pastoris DSMZ 70382",
		"Kop_b" => "Pichia pastoris GS115",
		"Lak_b" => "Saccharomyces kluyveri CBS 3082",
		"Lak_a" => "Saccharomyces kluyveri NRRL Y-12651",
		"Lat" => "Kluyveromyces thermotolerans CBS 6340",
		"Lw" => "Kluyveromyces waltii NCYC 2644",
		"Low" => "Lodderomyces elongisporus NRLL YB-4239",
		"Mef" => "Metschnikowia fructicola ARO 277",
		"Mrg" => "Yamadazyma guilliermondii, Candida guilliermondii var. carpophila, Pichia guilliermondii ATCC 6260",
		"Mif" => "Pichia farinosa, Saccharomyces farinosus, Pichia farinosa var. farinosa, Yamadazyma farinosa, Candida cacaoi, Pichia sorbitophila",
		"Nac" => "Saccharomyces castellii Capriotti, Saccharomyces castellii, Naumovia castellii NRRL Y-12630",
		"Oga" => "Pichia angusta NCYC 495 leu1.1, Pichia angusta ATCC MYA-335",
		"Ogp" => "Hansenula polymorpha DL-1, Ogataea angusta DL-1 (misnomer), Ogataea parapolymorpha ATCC 26012, Pichia angusta DL-1",
		"Oap" => "Ogataea angusta CBS 4732, Pichia angusta CBS 4732",
		"Piu" => "Issatchenkia orientalis, Candida krusei, Saccharomyces sp. P1",
		"Phm" => "Pichia membranifaciens (E.C.Hansen) E.C.Hansen, Pichia membranaefaciens",
		"Sc_an" => "Hormiscium cerevisiae (obsolete), Mycotorula robusta (obsolete), Saccharomyces boulardii (obsolete), Saccharomyces carlsbergensis (obsolete), Saccharomyces oviformis, Saccharomyces italicus, Saccharomyces capensis, Saccharomyces uvarum var. melibiosus",
		"Sc_aq" => "Hormiscium cerevisiae (obsolete), Mycotorula robusta (obsolete), Saccharomyces boulardii (obsolete), Saccharomyces carlsbergensis (obsolete), Saccharomyces oviformis, Saccharomyces italicus, Saccharomyces capensis, Saccharomyces uvarum var. melibiosus",
		"Sc_al" => "Hormiscium cerevisiae (obsolete), Mycotorula robusta (obsolete), Saccharomyces boulardii (obsolete), Saccharomyces carlsbergensis (obsolete), Saccharomyces oviformis, Saccharomyces italicus, Saccharomyces capensis, Saccharomyces uvarum var. melibiosus",
		"Sc_b" => "Saccharomyces oviformis, Saccharomyces italicus, Saccharomyces capensis, Saccharomyces uvarum var. melibiosus",
		"Sc_c" => "Hormiscium cerevisiae (obsolete), Mycotorula robusta (obsolete), Saccharomyces boulardii (obsolete), Saccharomyces carlsbergensis (obsolete), Saccharomyces oviformis, Saccharomyces italicus, Saccharomyces capensis, Saccharomyces uvarum var. melibiosus",
		"Sc_a" => "Hormiscium cerevisiae (obsolete), Mycotorula robusta (obsolete), Saccharomyces boulardii (obsolete), Saccharomyces carlsbergensis (obsolete), Saccharomyces oviformis, Saccharomyces italicus, Saccharomyces capensis, Saccharomyces uvarum var. melibiosus",
		"Sc_am" => "Hormiscium cerevisiae (obsolete), Mycotorula robusta (obsolete), Saccharomyces boulardii (obsolete), Saccharomyces carlsbergensis (obsolete), Saccharomyces oviformis, Saccharomyces italicus, Saccharomyces capensis, Saccharomyces uvarum var. melibiosus",
		"Shc" => "Pichia stipitis CBS 6054",
		"Ttp" => "Kluyveromyces phaffii",
		"Tod" => "Saccharomyces fermentati, Saccharomyces rosei",
		"Vp" => "Vanderwaltozyma polyspora DSMZ 70294",
		"Wa" => "Pichia anomala NRRL Y-366",
		"Zr" => "Saccharomyces rouxii (obsolete)"
	}

	def species_table
		coldata = refdata2species
		content_tag(:thead, escape: false) do
			content_tag(:tr) do 
				content_tag(:th, "Species (teleomorph") +
				content_tag(:th, "Anamorph") +
				content_tag(:th, "Alternative names") +
				content_tag(:th, "Abbreviation")
			end
		end +
		content_tag(:tbody) do
			coldata.each.collect do |abbr, data|
				content_tag(:tr) do
					col1_safe = data[0]
					col2_safe = data[1]
					col3_safe = data[2]
					col4_safe = abbr.html_safe? ? abbr : abbr.html_safe
					content_tag(:td, col1_safe) + 
					content_tag(:td, col2_safe) + 
					content_tag(:td, col3_safe) +
					content_tag(:td, col4_safe)
				end
			end.join.html_safe
		end

	end


	def refdata2species
		ref_data = Helper.load_ref_data

		spabbr2xxx = {}

		ref_data.each do |key, data|
			genes = data["genes"]

			genes.each do |abbr, genedat|
				sp_abbr = abbr.match(/
				([A-Z][_a-z]+) # match species abbreviation -- everything until the 2. uppercase char
				[A-Z][\w]+ # protein class and variant -- everything after 2. uppercase char 
				/x)[1]

				if ! spabbr2xxx.has_key?(sp_abbr) then
					begin
						anamorph = Anamorphs[sp_abbr] || ""
						altname = Altnames[sp_abbr] || ""
						name = genedat["species"]
					rescue
						anamorph = "?" if ! anamorph
						altname = "?" if ! altname
						name = "?" if ! name
					ensure
					spabbr2xxx[sp_abbr] = [name.html_safe, anamorph.html_safe, altname.html_safe]
					end
				end
			end
		end

		return spabbr2xxx

	end
end
