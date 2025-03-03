namespace Sequence.Position.Extract
{
    /// <summary>
    /// 同じ染色体位置を持つ項目のコレクション
    /// </summary>
    /// <typeparam name="T">ゲノム位置情報を持つクラス</typeparam>
    internal sealed class SameChrGenomePositionItems<T>
        where T : IHasGenomePositionItem
    {
        private static readonly GenomePositionComparer _genomePositionComparer = new();

        private readonly T[] _items;
        private readonly int _chunkSize;

        /// <summary>
        /// 同じ染色体位置を持つ項目のコレクションを作成する。
        /// </summary>
        /// <param name="values">同じ染色体位置を持つ項目</param>
        /// <param name="chunkSize">chunkの大きさ(0以下なら自動調整)</param>
        public SameChrGenomePositionItems(IEnumerable<T> values, int chunkSize)
        {
            _items = [.. values];
            _chunkSize = chunkSize;

            var uniqChrCount = _items.Select(x => x.GenomePosition.ChrName).Distinct().Count();
            if (uniqChrCount != 1) throw new ArgumentException("異なる染色体があります。");
        }

        /// <summary>
        /// 所持している項目を一定個数毎に分ける。
        /// </summary>
        /// <returns>ゲノム位置項目の配列</returns>
        public GenomePositionChunk<T>[] Chunk()
        {
            var chunkSize = GetChunkSize();

            var posChunks = _items
                .OrderBy(x => x.GenomePosition, _genomePositionComparer)
                .Chunk(chunkSize)
                .Select(x => new GenomePositionChunk<T>(x))
                .ToArray();

            return posChunks;
        }

        /// <summary>
        /// chunkサイズを取得する。
        /// _chunkSizeが0以下の場合は項目数をもとに算出する。
        /// </summary>
        /// <returns>chunkSize</returns>
        private int GetChunkSize()
        {
            if(_chunkSize > 0) return _chunkSize;

            return (int)Math.Sqrt(_items.Length) + 1;
        }
    }
}
