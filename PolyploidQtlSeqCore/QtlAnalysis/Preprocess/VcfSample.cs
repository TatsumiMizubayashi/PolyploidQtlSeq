using PolyploidQtlSeqCore.QtlAnalysis.Preprocess.IO;

namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess
{
    /// <summary>
    /// VCFに記載されているサンプル情報
    /// </summary>
    internal class VcfSample
    {
        /// <summary>
        /// VCFに記載されているサンプル情報を作成する。
        /// </summary>
        /// <param name="gt">GT</param>
        /// <param name="ad">AD</param>
        /// <param name="allele">アレル</param>
        public VcfSample(GT gt, AD ad, string allele)
        {
            GtType = gt.Type;
            Allele = allele;
            RefCount = ad.RefCount;
            AltCount = ad.AltCount;
            Depth = ad.Depth;
        }

        /// <summary>
        /// GTの種類を取得する。
        /// </summary>
        public GtType GtType { get; }

        /// <summary>
        /// アレルを取得する。
        /// </summary>
        public string Allele { get; }

        /// <summary>
        /// Depthを取得する。
        /// </summary>
        public int Depth { get; }

        /// <summary>
        /// REF型リード数を取得する。
        /// </summary>
        public int RefCount { get; }

        /// <summary>
        /// ALT型リード数を取得する。
        /// </summary>
        public int AltCount { get; }
    }
}
