namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess.IO
{
    /// <summary>
    /// Multiアレルタイプのリードカウンター
    /// </summary>
    internal class MultiAlleleReadCounter : ReadCounter
    {
        /*
         * Multiの場合は解析対象外なので、後に除かれる。
         * 何らかの値を返す必要があるので以下のルールを適応する。
         * refCountはRef型リード数
         * altCountはAlt型で最大リード
         */


        /// <summary>
        /// リード数をカウントする。
        /// </summary>
        /// <param name="ad">AD</param>
        /// <returns>(refCount, altCount)</returns>
        public override (int refCount, int altCount) Count(AD ad)
        {
            var refCount = ad.ReadCounts[0];
            var altCount = ad.ReadCounts.Skip(1).Max();

            return (refCount, altCount);
        }
    }
}
