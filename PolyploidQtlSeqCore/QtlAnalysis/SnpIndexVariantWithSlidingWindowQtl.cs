using PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow;
using Sequence.Position;
using Sequence.Position.Extract;

namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// SlidingWindowのQTL情報を持ったSnpIndex変異
    /// </summary>
    internal class SnpIndexVariantWithSlidingWindowQtl : IHasGenomePositionItem
    {
        /// <summary>
        /// SlidingWindowのQTL情報を持ったSnpIndex変異を作成する。
        /// </summary>
        /// <param name="variant">QTL情報を持つ変異情報</param>
        /// <param name="maxScoreWindowQtl">最大スコアWindowQTL情報</param>
        public SnpIndexVariantWithSlidingWindowQtl(SnpIndexVariantWithQtl variant, MaxScoreSlidingWindowQtl maxScoreWindowQtl)
        {
            GenomePosition = variant.GenomePosition;
            RefAllele = variant.RefAllele;
            Type = variant.Type;
            Annotations = variant.Annotations;
            Parent1 = variant.Parent1;
            Bulk1 = variant.Bulk1;
            Bulk2 = variant.Bulk2;
            DeltaSnpIndex = variant.DeltaSnpIndex;
            P95ThresholdDeltaSnpIndex = variant.P95ThresholdDeltaSnpIndex;
            P99ThresholdDeltaSnpIndex = variant.P99ThresholdDeltaSnpIndex;
            PValue = variant.PValue;
            P95Qtl = variant.P95Qtl;
            P99Qtl = variant.P99Qtl;
            Score = variant.Score;
            MaxScoreSlidingWindowQtl = maxScoreWindowQtl;
        }

        #region SnpIndexVariantWithQtlからの持ち込み

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

        #endregion

        /// <summary>
        /// 最大スコアWindowQTL情報を取得する。
        /// </summary>
        public MaxScoreSlidingWindowQtl MaxScoreSlidingWindowQtl { get; }

    }
}
