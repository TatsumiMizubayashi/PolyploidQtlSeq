using Sequence.Position;
using Sequence.Position.Extract;

namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess
{
    /*
     * VCFファイルから取得した必要な情報を入れるだけの器。
     */

    /// <summary>
    /// VCFに記載されている変異情報
    /// </summary>
    internal class VcfVariant : IHasGenomePositionItem
    {
        /// <summary>
        /// VCF変異情報を作成する。
        /// </summary>
        /// <param name="genomePosition">ゲノム位置</param>
        /// <param name="refAllele">Refアレル塩基</param>
        /// <param name="type">変異の種類</param>
        /// <param name="isMultiAltAllele">ALTが複数あるかどうか</param>
        /// <param name="parent1">親1</param>
        /// <param name="parent2">親2</param>
        /// <param name="bulk1">Bulk1</param>
        /// <param name="bulk2">Bulk2</param>
        /// <param name="annotations">SnpEffアノテーション</param>
        public VcfVariant(GenomePosition genomePosition, string refAllele, VariantType type, bool isMultiAltAllele,
            VcfParent1 parent1, VcfParent2 parent2, VcfBulk1 bulk1, VcfBulk2 bulk2, SnpEffAnnotations annotations)
        {
            GenomePosition = genomePosition;
            RefAllele = refAllele;
            Type = type;
            IsMultiAltAllele = isMultiAltAllele;
            Parent1 = parent1;
            Parent2 = parent2;
            Bulk1 = bulk1;
            Bulk2 = bulk2;
            Annotations = annotations;
        }

        /// <summary>
        /// ゲノム位置を取得する。
        /// </summary>
        public GenomePosition GenomePosition { get; }

        /// <summary>
        /// Refアレルの塩基を取得する。
        /// </summary>
        public string RefAllele { get; }

        /// <summary>
        /// 変異の種類を取得する。
        /// </summary>
        public VariantType Type { get; }

        /// <summary>
        /// Altアレルが複数あるかどうかを取得する。
        /// </summary>
        public bool IsMultiAltAllele { get; }

        /// <summary>
        /// 親1の情報を取得する。
        /// </summary>
        public VcfParent1 Parent1 { get; }

        /// <summary>
        /// 親2の情報を取得する。
        /// </summary>
        public VcfParent2 Parent2 { get; }

        /// <summary>
        /// Bulk1の情報を取得する。
        /// </summary>
        public VcfBulk1 Bulk1 { get; }

        /// <summary>
        /// Bulk2の情報を取得する。
        /// </summary>
        public VcfBulk2 Bulk2 { get; }

        /// <summary>
        /// SnpEffアノテーションコレクションを取得する。
        /// </summary>
        public SnpEffAnnotations Annotations { get; }
    }
}
