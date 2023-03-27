using PolyploidQtlSeqCore.QtlAnalysis.Preprocess;
using Sequence.Position;
using Sequence.Position.Extract;

namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// SNP-indexを持つ変異情報
    /// </summary>
    internal class SnpIndexVariant : IHasGenomePositionItem
    {
        /// <summary>
        /// SNP-indexを持つ変異情報を作成する。
        /// </summary>
        /// <param name="vcfVariant"></param>
        public SnpIndexVariant(VcfVariant vcfVariant)
        {
            GenomePosition = vcfVariant.GenomePosition;
            RefAllele = vcfVariant.RefAllele;
            Type = vcfVariant.Type;
            Annotations = vcfVariant.Annotations;

            Parent1 = new Parent1(vcfVariant.Parent1);
            Parent2 = new Parent2(vcfVariant.Parent1, vcfVariant.Parent2);
            Bulk1 = new Bulk1(vcfVariant.Parent1, vcfVariant.Bulk1);
            Bulk2 = new Bulk2(vcfVariant.Parent1, vcfVariant.Bulk2);
            DeltaSnpIndex = Bulk2.CalcDeltaSnpIndex(Bulk1);
        }

        /// <summary>
        /// ゲノム位置を取得する。
        /// </summary>
        public GenomePosition GenomePosition { get; }

        /// <summary>
        /// Refアレル塩基を取得する。
        /// </summary>
        public string RefAllele { get; }

        /// <summary>
        /// 変異の種類を取得する。
        /// </summary>
        public VariantType Type { get; }

        /// <summary>
        /// SnpEffアノテーションコレクションを取得する。
        /// </summary>
        public SnpEffAnnotations Annotations { get; }

        /// <summary>
        /// 親1情報を取得する。
        /// </summary>
        public Parent1 Parent1 { get; }

        /// <summary>
        /// 親2情報を取得する。
        /// </summary>
        public Parent2 Parent2 { get; }

        /// <summary>
        /// Bulk1情報を取得する。
        /// </summary>
        public Bulk1 Bulk1 { get; }

        /// <summary>
        /// Bulk2情報を取得する。
        /// </summary>
        public Bulk2 Bulk2 { get; }

        /// <summary>
        /// ΔSNP-indexを取得する。
        /// </summary>
        public DeltaSnpIndex DeltaSnpIndex { get; }
    }
}
