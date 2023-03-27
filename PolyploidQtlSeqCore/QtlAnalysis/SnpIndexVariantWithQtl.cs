using Sequence.Position;
using Sequence.Position.Extract;

namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// 変異のQTL解析情報を持つSnpIndex変異情報
    /// </summary>
    internal class SnpIndexVariantWithQtl : IHasGenomePositionItem
    {
        /// <summary>
        /// 変異のQTL解析情報を持つSnpIndex変異情報を作成する。
        /// </summary>
        /// <param name="variant">SNP-index変異</param>
        /// <param name="p95Threshold">P95しきい値</param>
        /// <param name="p99Threshold">P99しきい値</param>
        /// <param name="pValue">P値</param>
        public SnpIndexVariantWithQtl(SnpIndexVariant variant, QtlThresholdDeltaSnpIndex p95Threshold,
            QtlThresholdDeltaSnpIndex p99Threshold, PValue pValue)
        {
            GenomePosition = variant.GenomePosition;
            RefAllele = variant.RefAllele;
            Type = variant.Type;
            Annotations = variant.Annotations;
            Parent1 = variant.Parent1;
            Bulk1 = variant.Bulk1;
            Bulk2 = variant.Bulk2;
            DeltaSnpIndex = variant.DeltaSnpIndex;

            P95ThresholdDeltaSnpIndex = p95Threshold;
            P99ThresholdDeltaSnpIndex = p99Threshold;
            PValue = pValue;

            P95Qtl = DeltaSnpIndex.IsQtl(P95ThresholdDeltaSnpIndex);
            P99Qtl = DeltaSnpIndex.IsQtl(P99ThresholdDeltaSnpIndex);
            Score = PValue.ToScore();
        }

        #region SnpIndexVariantから持ち込まれた情報
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

        #endregion

        /// <summary>
        /// P95のΔSNP-indexしきい値を取得する。
        /// </summary>
        public QtlThresholdDeltaSnpIndex P95ThresholdDeltaSnpIndex { get; }

        /// <summary>
        /// P99のΔSNP-indexしきい値を取得する。
        /// </summary>
        public QtlThresholdDeltaSnpIndex P99ThresholdDeltaSnpIndex { get; }

        /// <summary>
        /// P値を取得する。
        /// </summary>
        public PValue PValue { get; }

        /// <summary>
        /// P95でQTLがあるかどうかを取得する。
        /// </summary>
        public bool P95Qtl { get; }

        /// <summary>
        /// P99でQTLがあるかどうかを取得する。
        /// </summary>
        public bool P99Qtl { get; }

        /// <summary>
        /// スコア(-log10(P))を取得する。
        /// </summary>
        public Score Score { get; }
    }
}
