namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess.IO
{
    /// <summary>
    /// RawVCFサンプル
    /// </summary>
    internal class RawVcfSample
    {
        /// <summary>
        /// VcfSampleインスタンスを作成する。
        /// </summary>
        /// <param name="ad">AD</param>
        /// <param name="gt">GT</param>
        public RawVcfSample(AD ad, GT gt)
        {
            AD = ad;
            GT = gt;
        }

        /// <summary>
        /// ADを取得する。
        /// </summary>
        public AD AD { get; }

        /// <summary>
        /// GTを取得する。
        /// </summary>
        public GT GT { get; }

        /// <summary>
        /// ADを更新したVcfSampleを作成する。
        /// </summary>
        /// <param name="readCounter">ReadCounter</param>
        /// <returns>更新したVcfSample</returns>
        public RawVcfSample Update(ReadCounter readCounter)
        {
            var updateAd = new AD(AD, readCounter);

            return new RawVcfSample(updateAd, GT);
        }
    }
}
