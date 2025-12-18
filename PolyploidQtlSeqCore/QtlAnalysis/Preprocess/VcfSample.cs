using PolyploidQtlSeqCore.QtlAnalysis.Preprocess.IO;

namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess
{
    /// <summary>
    /// VCFに記載されているサンプル情報
    /// </summary>
    internal class VcfSample
    {
        private static readonly char[] _splitter = ['/', '|'];
        private const string _noData = ".";
        private const string _refValue = "0";

        /// <summary>
        /// VCFに記載されているサンプル情報を作成する。
        /// </summary>
        /// <param name="gt">GT値</param>
        /// <param name="allele">アレル</param>
        /// <param name="refCount">Refリード数</param>
        /// <param name="altCount">Altリード数</param>
        [Obsolete]
        public VcfSample(string gt, string allele, int refCount, int altCount)
        {
            GtType = ToGtType(gt);
            Allele = allele;
            Depth = refCount + altCount;
            RefCount = refCount;
            AltCount = altCount;
        }

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

        private static GtType ToGtType(string value)
        {
            var items = value.Split(_splitter).Distinct().ToArray();
            if (items.Length == 2) return GtType.Hetero;

            var gtValue = items[0];
            if (gtValue == _noData) return GtType.NoData;
            if (gtValue == _refValue) return GtType.RefHomo;

            return GtType.AltHomo;
        }
    }
}
